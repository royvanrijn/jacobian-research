# The abc exceptional set and powerful triples

For \(\lambda>0\), let \(N_\lambda(X)\) count ordered positive primitive triples \(a+b=c\le X\) with \(\operatorname{rad}(abc)<c^\lambda\). The \(abc\) conjecture is the assertion that \(N_\lambda(X)\) is bounded in \(X\) for every \(0<\lambda<1\). Let \(T_{p,q,r}(X)\) count ordered positive primitive triples \(a+b=c\le X\) for which the three terms are respectively \(p\)-full, \(q\)-full and \(r\)-full. An integer is \(s\)-full when every positive prime valuation is at least \(s\), so \(1\) is \(s\)-full for every \(s\). [Bernert–Browning–Lichtman–Teräväinen, §1](https://arxiv.org/html/2410.12234v2#S1) · [Browning–Verzobio, §1](https://arxiv.org/html/2608.24512v1#S1)

The repository's unrestricted small-budget estimate has limiting
coefficient \(1/3\). At unit budget the later two-tall-term
certificates give exponent \(1070000267/2000000686\).
See the [small-budget theorem](ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13),
the [fixed-budget bounds](#exceptional-set-bounds-from-the-two-tall-terms),
and the [strict high-powerful transfer](POWERFUL_TRIPLES_ABC_ENDPOINT_TRANSFER.md#theorem-a-strict-saving-from-endpoint-diffusion)
for their separate ranges.

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

Hanxin Zhang's August preprint gives the claimed uniform
exponent \(34/59\) for \(0<\lambda\le1\), retaining all
valuation layers. Its exact all-depth checker replays at
the pinned source revision. The new prime-power fibre
estimate adds an arithmetic row on divisibility sums
rather than fixed small valuation indices. It excludes
mass concentrated at a common high prime degree even
when that degree lies outside every fixed low-index
selector. Distinct prime-degree profiles still evade this
row family.
[All-depth certificate and fibre transfer](ABC_ALL_DEPTH_VALUATION_TRANSFER.md)

Grouping divisible layers extends the exposed fourth-moment
estimate to a sparse root dictionary of larger height.
The new mixed \(L^2,L^4,L^4\) rows strengthen several
concentrated valuation profiles. The audited external
rows, even with moving selectors, still admit
formal narrow prime-degree profiles near \(2\lambda/3\).
In the same-target collision argument the diagonal arm
alone is at least \(\lambda/3\) on balanced data.
[Grouped moments](ABC_GROUPED_VALUATION_FOURTH_MOMENTS.md)
· [Moving-selector limitation](ABC_EXTERNAL_MOVING_ROWS_LIMIT.md)
· [Current source versions, closed-budget BBLT deduction and Pasten inversion](ABC_PRIMARY_BOUNDS_AND_POINTWISE_BRIDGES.md)


Pasten's conductor-dependent Faltings height bound is
\(h(E)\ll M\phi(D)\log\ell\) for an admissible
factorization \(N_E=DM\) and a prime \(\ell\nmid N_E\).
It gives \(h(E)\ll N_E\log\log N_E\) in general and
\(h(E)\ll_S N_E\) for curves semistable outside a fixed
set \(S\). These are bounds on height at a given conductor.
The introduction records the stronger \(N_E^{1/3+\varepsilon}\)
height bound in the special Frey setting. Our relevant
quantity here is the number of occurring support records
at polynomial radical budget. Bilu and Luca prove finite
exceptions to the common-torsion explanation for
multiplicatively dependent pairs of nonzero sums from two
fixed, almost disjoint finitely generated multiplicative
groups. Their group hypotheses are fixed, whereas the
supports in our secant and transfer dictionaries can vary
through the height limit.
[Pasten, Theorem 1.1 and Corollaries 1.2–1.3](https://arxiv.org/html/2609.17390v1#S1)
· [Bilu–Luca, Theorem 1.2](https://arxiv.org/html/2607.28857v1#S1.SS1)

Chan and Koymans count integral points near the square–cube cusp by
using continued fractions in a transformed quartic. Their
Theorem 1.2 is unconditional and counts by discriminant size. For
the primitive split Frey model \(y^2=x(x-a)(x+b)\), its invariants
before minimalization satisfy
\[
 c_4=16(a^2+ab+b^2),\qquad
 \Delta=16(abc)^2,\qquad
 \frac{1728|\Delta|}{c_4^3}\asymp
        \left(\frac{\min(a,b)}{c}\right)^2 .
\]
The ratio is invariant under minimalization. At the minimum
nonatomic endpoint all three original numerical height exponents
tend to one. Thus this ratio can have exponent zero, and a
small radical of the discriminant alone does not place the Frey
invariants in a power-thin numerical cusp. This identifies the
height input missing from a direct transfer of the cusp estimate.
[Chan–Koymans, Theorem 1.2 and §§2–3](https://arxiv.org/html/2611.31390v1#S1)

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

The small-determinant saving also holds for arbitrary nonzero
determinants of absolute size at most \(c^{\eta_0\lambda^3}\),
with the same type of bounds on the two multipliers and on
the product of odd overlaps. A uniform estimate for coprime
pairs at small distance extends Diao's balanced determinant
argument and the unbalanced quadratic-field argument.
When the quadratic norm splits, factoring the difference
of squares preserves the squarefree kernel product and
halves the numerical height exponent, so the split branch
descends. The resulting small-budget saving is an absolute
multiple of \(\lambda^3\) on this sector.
[Small nonunit determinants](ABC_SMALL_NONUNIT_DETERMINANTS.md)

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

The one-third argument also counts by an effective radical
budget. In the diagonal range of the residue moments, the
modulus threshold uses actual root sizes and the count uses
root dictionary cardinalities. It follows that a sum of
dictionary exponents \(B\) gives the bound
\(X^{B/3+o(\lambda)+\varepsilon}\) in the iterated small-budget
limit. Prime prices with a sublogarithmic partition sum
convert this to a global count. Any fixed positive fraction
of radical weight on a prescribed pool of polylogarithmic
cardinality suffices for a strict saving, with no fixed
upper limit on the logarithmic power of its cardinality.
The weighted two-letter inequality gives a fixed-value
fibre exponent \((1/6+o(1))\lambda\) at balance. In general
the pointwise exponent is at most one quarter of the
radical budget outside the prescribed primitive letter,
up to the usual small-budget allowance. This also bounds
the pointwise count at a prescribed least-repaired target.
[Fixed complete values and prescribed least-repaired targets](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md#fixed-complete-values-and-prescribed-least-repaired-targets)

At saturation, the logarithmic radical weight is
asymptotically concentrated on primes whose occurrence
frequencies are \(p^{-1+o(1)}\) in the weighted mean.
Anatomic root marginals attain the sum of their numerical
height caps in entropy. Write \(W_T\) for the charged
joint dictionary exponent on letter \(T\). In original
height shells the improved count is
\[
                      D\le \tfrac12\bigl(W-\max_T W_T\bigr)+o(\lambda),
                       \qquad W=\sum_T W_T .
\]
The allowed within-letter groups have bounded size or
boundedly many shrinking valuation windows. Already
frozen roots retain their raw caps and have zero
conditional dictionary cost. Equality at \(W/3\) forces
balance of all three effective costs. A zero numerical
short height is handled by releasing the tall pair in the
real short strip. [Effective radical budget](ABC_EFFECTIVE_RADICAL_BUDGET.md)
· [Grouped root dictionaries](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md)

A prescribed cheap full lock also permits direct signed
ratio moments on dictionaries spanning the two free
letters. For any partition into bounded root groups,
their total joint entropy saturates the two free-letter
cap sum. Groups spread across a bounded number of
shrinking valuation windows may be arbitrarily large at
fixed truncation, by the wide-cofactor count. Two
simultaneous locks permit dictionaries spanning all three
letters and bound their joint cap sum from below by
\(2\lambda/3\) at saturation. A fixed real height gap gives the tall-pair bound without a
lock and permits three-letter groups with only one full lock
on a tall letter. Abstract modular dictionaries can attain the
square-root count while satisfying the signed diagonal moment
bounds on arbitrary unions of groups. An arithmetic improvement
has to use restrictions on the actual values beyond those
moment bounds and the local pair entropies.
[Real and modular joint dictionaries](ABC_REAL_STRIP_GROUP_DICTIONARIES.md)
· [Signed residue moment obstruction](ABC_THREE_RESIDUE_DICTIONARY_MOMENT_OBSTRUCTION.md)

The effective budget also accepts conditional prices
along ordered prime supports. A sublogarithmic next-symbol
partition sum holds uniformly at every prefix. This
turns successive conditional prime frequencies into
admissible dictionary costs. At saturation their
logarithmic radical-weighted deviations from \(1/p\)
again vanish in mean, within every block covered by the
relevant moment bound. A quantitative application
saves when a fixed positive share of the budget is
carried by close-prime pairs within those blocks,
without one global prescribed prime pool.
[Across-letter dictionaries and simultaneous moduli](ABC_CROSS_LETTER_ROOT_DICTIONARY_ENTROPY.md)
· [Conditional prime-support coding](ABC_CONDITIONAL_PRIME_SUPPORT_CODING.md)


An actual two-letter power norm can count distinct targets at
equal heights. Fix full divisors \(F,G\) of the two summands and
retain a target packet \(T=C R^k\). For a degree-\(q\) binomial norm
the finite-place embedding determinant uses precision
\((FG)^{(q-1)/2}C_0^{\binom{q-1}{2}}\), where \(C_0\) is the
part of \(C\) at its good primes. The archimedean size is controlled
by \(R\). Already paid coefficient supports can contribute their
full locks. The resulting spherical \((2,2,4)\) sector has
exponent \((1/3-1/100)\lambda+\varepsilon\) under the stated
caps. Dense prime valuation indices leave the two binomial
coefficient records too expensive.
[Finite-place norm target counts](ABC_FINITE_PLACE_NORM_TARGET_COUNTS.md)

The small-budget bound transfers in the reverse direction to
one-letter powered-number spacing. Normalizing a nearby pair in
a unit progression block charges only the radical of its step
in the derived triple. This gives quantitative local occupancy
and a global nearby-pair bound, with improved scaling recovery
for full integers. A second elementary argument stops a smooth
full divisor just above the block width. It gives subpower
occupancy on blocks of normalized length \(n^h\) from any
positive diffuse high-valuation component of the actual target,
even with macroscopic primes elsewhere. Both estimates still
permit one target in each of many different blocks.
[Powered progression occupancy](ABC_POWERED_PROGRESSION_OCCUPANCY.md)

An invertible \(p\)-adic Jacobian can be charged once by using
nearby zeros over the algebraic closure. For a square system,
the sum of its two largest Smith valuations is a sufficient
residual-depth threshold. A generic finite reference fibre
preserves the original congruences and allows residue labels.
In the primitive-kernel one-root argument this removes the
separate terminal inequality. The strengthened rank conditions
still fail on the symmetric diffuse test profile, across
arbitrary retained submeasures and fixed mesh choices.
[Multivariate approximate zeros](ABC_MULTIVARIATE_APPROXIMATE_ZEROES.md)
· [One-root target count and symmetric obstruction](ABC_ALGEBRAIC_ROOT_TARGET_OCCUPANCY.md)

Prescribing the actual parity defects on all three primitive letters
allows a second use of the affordable divisor records. The resulting
conic has a rational model of controlled height and each full power
gives half its numerical weight as parameter precision. The real
interval near a short letter supplies the matching factor of one half
in the determinant bound. This proves an absolute gap below \(1/3\)
for an absolute positive parity-defect budget, including arbitrary
nonatomic distributions of the even valuation degrees and macroscopic
primes. After the squarefree kernels have been counted, their support
gives the full odd-valuation weight at subpower cost.
[Square kernels and actual full divisors](ABC_SQUAREFREE_KERNEL_PROJECTIVE_LOCKS.md)

For fixed abelian varieties over varying bounded-degree prime-power
extensions, a canonical-height ball can be covered by smaller balls
whose differences lie in proper subfields. The discriminant controls
the spacing, and genus-theory rank bounds use each subfield's own
discriminant. This gives a subpower bound uniformly in the field,
and a count \(X^{f\lambda+\varepsilon}\) on the Euclidean signatures
\((3,3,3)\) and \((2,4,4)\) with joint defect radical at most
\(X^{f\lambda}\). Both conclusions apply directly to primitive
quotients under minimal repair. The fixed-base literature also gives
these two signature counts. A moving Legendre curve gives a separate
saving when a fourth determinant has a cheap squarefree kernel.
[Euclidean signatures and bounded-height points](ABC_EUCLIDEAN_GRID_POINT_COUNTS.md)

For varying shifts, two primitive columns in the same row with a
common target power record force its cofactor to divide one nonzero
difference determinant. A weighted power difference then gives
subpower column-pair codegree. The complete-label application raises
the required cofactor-dictionary exponent from \(1/12\) to \(1/6\)
at the minimum nonatomic endpoint. This estimate keeps the third
letter's arithmetic even when the shift is tall.
[Codegrees for varying shifts](ABC_SHIFT_DICTIONARY_WEIGHTED_CODEGREES.md)

High-order monomial coordinates in a nondegenerate relation can
be released even if they occur on several terms. A low-degree
component with a point at which exactly one coordinate vanishes
would have a function-field unit ratio of excessive order at
that point. When it has no such point, an arithmetic
Nullstellensatz fixes the prime support of each coordinate.
This gives a new mixed off-diagonal bound through prime-incidence
splitting and a fixed-target bound for powerful integers on a few
prescribed valuation layers. On the original equation it also
saves when more than two thirds of the radical budget is captured
by an arbitrary slowly growing set of high valuation positions. At the balanced endpoint the
diagonal still permits distinct targets with one pair each.
[Nondegenerate power relations and collisions](ABC_NONDEGENERATE_POWER_COLLISIONS.md)

The uniform missing estimate is still a power saving for the
distinct height-admissible targets of small radical in an affordable
two-divisor record. If its normalized cost is \(t\), a sufficient
conditional fibre bound is \(n^{1/3-t-\sigma}\) for some absolute
\(\sigma>0\) on one affordable pair in each balanced sequence,
with allowances smaller than \(\sigma\) fixed before the height
limit. The local graph and residue constraints alone leave this
possibility open. In particular, broad formal laws supported on
prime valuation indices have no cheap parity record or fixed
Euclidean defect record, and a target cofactor dictionary can
retain its full radical cap. These are profile obstructions to
the new estimates and do not supply primitive triples at
saturation.
[Affordable full-divisor records](ABC_WIDE_COFACTOR_CONGRUENCES.md#affordable-whole-coordinate-pair-locks-at-every-balanced-endpoint)

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

At a balanced endpoint the curve terminal can also start with
no paid full divisor. If the good-model exceptional support takes
away all nonatomic mass, the supported parts that were already
recorded supply a fixed positive-weight full modulus for the
purely atomic terminal. If nonatomic mass survives, the original
diffuse branch argument applies. Separately, a rational
parametrization of the complete integer values with at least
two distinct poles admits only subpower many parameters at
each fixed budget, by a uniform ideal-divisor and unit count.
These are per-curve counts. The number of curves needed for
an economical cover of an arbitrary pair fibre remains open.
[Curve terminals without a paid lock](ABC_CURVE_TERMINALS_WITHOUT_PAID_LOCKS.md)

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


#### Exact arithmetic sectors and determinant limits

An affordable record also couples distinct targets without
common-neighbor stars. Short integral additive relations on
the targets lift to exact relations on both summands. For
a fixed relative affine shape, counts for the normalized
patterns determine the three common contents by a rank-two
linear relation. The tuple count pays only the joint
defect-prime support. When two letters have normalized
progression products with no defects above the length \(k\),
Pratt's long-progression theorem bounds their shapes and
forces the third full divisor to be \(O(k(\log k)^3)\).
This excludes such configurations with a sufficiently tall
fixed divisor and leaves the unrestricted shape and defect
dictionaries open.
[Additive patterns and common contents](ABC_AFFORDABLE_FIBRE_ADDITIVE_PATTERNS.md)

A different constraint works on arbitrary sparse matchings.
After fixing concentrated supports, crossing two residual
summands gives a new primitive additive relation and a
normalized determinant. If the initial divisor-record cost is
\(t\), saturation requires an almost-all determinant radical
exponent at least \(2/3-4t\) when that is positive, up to
arbitrarily small allowances. A fixed target divisor further
divides each normalized determinant. Small quotient determinants
can therefore exclude some sparse fibres. Large-radical
determinants remain possible under this estimate.
[Determinants of distinct targets](ABC_AFFORDABLE_FIBRE_DETERMINANT_COMPLEXITY.md)

The size and valuation distribution of those determinants are
also constrained after thinning. A large class cannot have a
projective modulus within a fixed numerical margin of its
determinant height. Combining full factors from several
anchor secants, and optionally from the actual three query
letters, makes a thick set of high-order secant packets
incompatible with that margin at affordable root cost.
In the residue after removing the common modulus, almost
all numerical weight lies at bounded scaled orders on
primes outside both endpoint triples. A separate curve
count with finitely many anchor directions forces a
positive defect mass at those orders for every exact
power grid. The remaining obstruction includes large
radicals with low absolute valuations.
[Secant valuation weights](ABC_AFFORDABLE_FIBRE_SECANT_VALUATION_WEIGHTS.md)
· [Projective peeling and combined packets](ABC_AFFORDABLE_FIBRE_PROJECTIVE_SECANT_PEELING.md)
· [Power defects at bounded scaled orders](ABC_SECANT_VALUATION_DEFECT_GRAPHS.md#power-defects-at-bounded-scaled-orders-on-fresh-primes)

Bounded tuples in a large affordable fibre have exclusive
prime packets on every pairwise determinant. Each packet
retains almost all of the terminal angular numerical scale
and at least the residual fibre exponent in radical mass,
up to arbitrarily small allowances. Universal rational
transfers of projective directions with no new prime source
are cross ratios or, for multiple outputs, one common
change of markings. Each moving secant used by those
transfers contributes its exclusive packet to the resulting
primitive output. More general rational transfers with
mixed polynomial factors pay the fibre exponent for each
such factor on disjoint supports. The count uses
coefficient-uniform binary norm fibres in a nonlinear slot
and point-line incidences when every slot is linear.
[Cross ratios and fresh polynomial prime costs](ABC_CROSS_RATIO_FRESH_PRIME_OBSTRUCTION.md)

In the remaining two-mixed-factor pattern, dyadic
regularization accounts for the free-letter marginals.
A bounded sequence of support refinements of the
actual eligible tuple graph makes the original
pure-letter values recoverable from the transferred
denominator at low cost. The resulting output
expansion excludes a dense family with output budget
at most \(3r\) above a fixed power of \(n\). For
smaller outputs, the affordable congruences control
one-slot fibres for every same-orientation two-slot
polynomial chart and for opposite-orientation
two-slot additive charts. These bounds eliminate
the remaining two-slot additive transfer from
this route. A direct use for the original fibres
still requires an arithmetic upper bound on
transferred radicals.
[Pure-term slices](ABC_AFFORDABLE_PURE_TERM_SLICES.md)
· [Denominator recovery](ABC_PURE_TERM_DENOMINATOR_RECOVERY.md)


A three-divisor determinant gives a strict global count
for radicals largely supported on a prescribed
polylogarithmic prime pool. The primes in the pool may
themselves be large. The pool can have cardinality
\(O((\log X)^{\kappa_*})\) for an absolute \(\kappa_*>2\),
with absolute positive allowances for the outside radical
and for the exponent saving. Uniform pooled dispersion
lowers the weight-cost threshold, and the last pool prime
can be recorded separately even when it is macroscopic.
The resulting gap transfers to minimally repaired
primitive quotients in the corresponding sector.
[Polylogarithmic prime pools and divisor packing](ABC_SMOOTH_SUPPORT_DIVISOR_PACKING.md#uniform-saving-beyond-the-logarithmic-power-two)

For fixed polynomial models of local minimal repair, a
logarithmic derivative bound rules out a third-layer
polynomial divisibility identity when the three original
polynomials have all factor multiplicities at least three.
Integer parameters that pass the third-layer test despite
a nonzero additive error force a fixed factor to have
square-full values away from a fixed pool. Square and
cube twist bounds give a parameter saving for factors of
degree up to four, and higher depth handles any fixed
factor degree. This identifies an obstruction for the
polynomial model of the local tests.
[Polynomial models of third-layer repair](ABC_MINIMAL_REPAIR_POLYNOMIAL_THIRD_LAYER.md)

Complete odds quotients at three or four actual targets also
give a higher-order common-neighbor bound. Fixed relative odds
put the targets into independent linear forms in two primitive
common factors. A small mixed valuation defect bounds the
number of completions of those forms. Three actual targets
suffice at any fixed order \(q\ge3\), using a cyclic twist,
the primitive cubic count, or a quartic elliptic map.
Four actual targets suffice for parity, by the elliptic count
with rational \(2\)-torsion. The original summand kernels are
unrestricted in both target-only results. For example, if
the target \(q\)-defect has normalized mass at most \(1/72\)
and total nonatomic mass is at most \(1/2+1/60\), the
three-target count excludes the limiting coefficient \(1/3\).
The fixed-odds argument uses actual witnesses and their
target valuations. Bounding the mixed-defect dictionaries
on general stars remains a separate arithmetic step.
[Three-target bounds for fixed valuation orders](ABC_THREE_TARGET_PRIME_ORDER_STARS.md)
· [Four-target parity and mixed kernels](ABC_MULTI_TARGET_MIXED_KERNEL_COUNTS.md)

A common multiplier chosen to make all three letters exact powers
has a pointwise cost at every microscopic degree shared by the
three valuation spectra. With \(Q\) the least common multiple of
the chosen orders, the sum of the normalized log repaired height
\(\log(mc)/Q\) and the log coefficient radical is at least
\(\sum_e\min_i\log\operatorname{rad}_e(a_i)\). The minimum of the
joint cost is given by primewise CRT, including when the orders
are chosen after observing the triple. For pairwise coprime orders
an unrestricted height benchmark of \(Y^{1/Q}\) already occurs
in distinct minimally repaired heights of primitive squarefree
quotients. The two facts rule out a height-only transfer on
microscopically matching spectra even after removing arbitrary
common multipliers. Weakly matching degree limits can still have
small repair cost, by the disjoint-grid example.
[Exact-power repair obstruction](ABC_EXACT_POWER_REPAIR_OBSTRUCTION.md)
· [Coefficient records and the effective orders](ABC_EXACT_POWER_REPAIR_COEFFICIENTS.md)

For general orders the common-repair support condition is the
primitive grid given by the three effective shared orders. These
orders also give the ramification indices of the normalized
weighted Fermat cover. In particular large pairwise coprime
orders give the projective line, while shared orders impose
actual degree divisibility away from the coefficient primes.
This calculation leaves the one-target minimal-fullness
equivalence (MRG) as a distinct problem. Its restricted count
uses the primitive summand radicals along with the repaired
target height.
[Weighted Fermat cover](ABC_EXACT_POWER_REPAIR_COEFFICIENTS.md#the-effective-orders-determine-the-weighted-fermat-curve)
· [One-target minimal repair](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)

For minimally repaired primitive quotients, actual cyclotomic
splitting gives a saving using all three radical bounds.
Opposite primitive \(q\)-th powers yield the coefficient
\(1/(\tau(q)+2)\), with \(1/4\) for squares and \(1/6\) for
sixth powers. The fixed-factor count transfers the norm
recovery used for Erdős #943. For same-side weighted powers
it also gives a direct repaired-height bound when the target
binomial splits into several nonlinear rational factors.
A Gaussian norm family with trivial squarefree kernels
attains the full repaired exponent when the original summand
radical cap is dropped, even with comparable heights and a
wholly diffuse target spectrum.
[Cyclotomic factors and minimal repair](ABC_MINIMAL_REPAIR_CYCLOTOMIC_SPLIT.md)

For genuine square kernels with one short primitive summand,
angular packing in the imaginary quadratic field saves on
a near-dyadic-power target. The dyadic order permits a
subpower bound on the eligible ideal classes, so the count
pays only the summand kernel and target coefficient radicals.
[Angular packing of primitive square kernels](ABC_MINIMAL_REPAIR_CYCLOTOMIC_SPLIT.md#angular-packing-of-primitive-square-kernels-and-near-power-targets)


Signed valuation phases allow savings without exact power
divisibility on either free letter. The two basis variables
can overlap in prime support and the count applies to distinct
exact rational outputs. This excludes phase concentration on
more than \(1/3\) of the two free masses when the remaining
letter has positive nonatomic mass, and covers logarithmically
many genuinely short degree windows. A microscopic entropy
bound explains why the same small-height-product condition
does not follow from the weak degree measures.
[Signed valuation phases](ABC_SIGNED_VALUATION_PHASES.md)

Several genuinely unrelated grids can also be released.
A short kernel vector in the actual valuation strip reduces
a bounded-rank signed progression to two bases. If the
product of digit widths is \(o(1/L)\), the resulting
two-coordinate height product is \(o(L)\) even for
mismatched steps and offsets. Sparsity of the degree
support alone is insufficient: polynomially sparse diffuse
profiles can avoid positive mass on every such bounded-rank
grid.
[Sparse valuation grid compression](ABC_SPARSE_VALUATION_GRID_COMPRESSION.md)

Independent weighted difference entropy gives a common
grid criterion even when the literal degree supports
are full. A bounded-information coupling of two draws
does not give the same conclusion, as shown by two
coprime-step grids meeting only the endpoint release
mass. On the complementary letter, selecting actual
small primes inside high valuation layers gives a cheap
full modulus even when its scaled degree law is
atomic. Together with two shrinking free-letter degree
bands this gives limiting normalized exponent zero
on that sector.
[Entropy and full divisors inside a layer](ABC_VALUATION_ENTROPY_GRID.md)

A direct anisotropic rational congruence count also
allows a growing number of signed bases on polynomially
sparse width products. It uses a triangular derivative
descent through a bounded number of equations at each
fixed budget, retaining the same actual full modulus.
[Anisotropic monomial congruences](ABC_ANISOTROPIC_RATIONAL_CONGRUENCES.md)



With independent exact valuation orders on just two letters,
the count has exponent \((f+r/2+\eta)\lambda+\varepsilon\).
Here \(f\lambda\) caps the joint coefficient radical on those
letters and \(r\lambda\) caps the radical of the third. The
orders may vary, with product at least a large fixed multiple
of \(1/\lambda\). A full modulus and an anisotropic determinant
handle the case without a large high-order prime. A single
macroscopic prime gives the third variable of a generalized
Fermat surface in the complementary case. At saturation
this bounds the sum of two actual divisibility masses by \(1/2\),
or \(1/3\) whenever the third spectrum contains a positive
nonatomic part. Using the bounds in all three orientations
excludes joint coefficient radical exponent below \(5\lambda/12\)
by a fixed margin when all three letters admit high exact orders.
For orders at least \(\lambda^{-\alpha}\), \(\alpha>2/3\), a thin-band
modulus raises the threshold to \(4\lambda/9\).
[Two valuation grids](ABC_TWO_VALUATION_GRIDS.md)


A small residual coefficient radical gives a \(7\lambda/24\)
exponent for triples with fine exact valuation divisibility on
all three letters. The grid has vanishing scaled mesh and so
allows wholly diffuse limiting profiles. On variable-shift
fibres, near-product density yields fixed multiplicative rows
whose full target divisors can be combined into a single CRT
modulus. Its numerical weight can exceed the actual shift width
even where a fixed-shift estimate gives no saving. The general
affordable-record count still needs a bound that uses the small
radical of the third letter across shifts. A separate affine bound
controls the occupancy of many short-interval fibres by their
smaller residual nonatomic mass, including sparse matchings.
[Affine count for sparse interval fibres](ABC_SPARSE_AFFINE_CONTENT.md)
· [Fine valuation divisibility](ABC_FINE_VALUATION_DIVISIBILITY.md)
· [Common rows with a varying shift](ABC_AGGREGATE_SHIFT_DICTIONARY_STARS.md)

A cyclic-cover descent gives a uniform count for two
weighted \(d\)-th powers and a third weighted
\(\ell\)-th power when \(d=\ell^j\) and the cover has
genus at least two. The point count in each primitive
coefficient fibre grows at most exponentially in the
number of coefficient primes. Summing coefficient defects
costs just their joint radical exponent \(f+h\). This
improves the repaired count whenever that exponent is
below the balanced target mass. It includes two exact
eighth-power grids with a square kernel on the third
letter, and does not require a height gap.
[Prime-power Fermat descent and divisibility masses](ABC_PRIME_POWER_FERMAT_DESCENT.md)

The three-point cover argument gives a height-free
polynomial bound in the defect radical for arbitrary fixed
hyperbolic signatures. Discriminant multiplicity would
make the support cost subpower by bounding both the number
of specialization fields and their class group torsion.
For a hyperbolic prime-primary part of the signature,
a cover with prime-group monodromy gives this stronger
bound unconditionally, including \((2,4,8)\), \((4,4,4)\)
and \((3,3,9)\). Pairwise coprime signatures have perfect
triangle group and admit no nontrivial solvable three-point
monodromy. The fixed-signature polynomial bound still
gives a subpower height count for primitive exact-power
Fermat–Catalan triples with all orders variable.
[Mixed valuation orders and specialization fields](ABC_MISMATCHED_GRID_NUMBER_FIELDS.md)

On dense consecutive valuation profiles, a regular
three-point cover pays at least \(4L/9-o(L)\) in
good-prime root discriminant, with equality in the
formal optimization at signature \((2,4,6)\).
For prime-degree profiles the lower charge is
\(85L/126-o(L)\). Changing inertia orders does not
remove this cost, and a moving twist producing an
unramified point must carry the same local inertia.
[Tame conductor bounds and the twisting cost](ABC_TAME_CONDUCTOR_PROFILES.md)

For one-target minimal repair, primitive row content can absorb
the entire extra common-slope Taylor divisor. Actual primitive
examples attain precisely the original Taylor valuation after
that content is removed. On the broad diffuse profile the exact
one-root interpolation test also excludes every support with
effective horizontal area at most one, allowing nonnested
supports. Paid subdivisions of only the narrow real variables
can be matched by finite-place depth on a diffuse modulus.
[Primitive content after repair](ABC_MINIMAL_REPAIR_PROJECTIVE_CONTENT.md)
· [Effective area for exact Taylor supports](ABC_EXACT_TAYLOR_EFFECTIVE_AREA.md)
· [The cost of real subdivision](ABC_WEIGHTED_HILBERT_TAYLOR.md#paying-for-short-real-intervals-by-additional-finite-place-precision)


Occupied-interval counts are not forced by overlapping
records, short difference determinants, or simultaneous record
entropy using only the two free radical caps. One primitive
matching attains the unsaved occupied-interval count for
every affordable slot record in those relaxed hypotheses.
Any global use of the interval bounds must retain the
small radical of the third letter across the distinct targets.
[Separated-interval obstruction](ABC_OVERLAPPING_FULL_RECORD_COUNTEREXAMPLES.md#simultaneous-separation-of-target-intervals)


For minimally repaired quotients, the genuine weighted
binomial norm can be counted in a narrow real chart even
when it is irreducible and all ideal classes vary with the
coefficients. In a pure prime-power field the matching
class torsion is subpower. The full embedding determinant
also gives a count after paying for real charts when neither
summand is short. The distinct quadratic norm split uses
a square norm defect and an actual reference solution; it
removes the reference ideal from the radical count and
does not require a short summand.
[Primitive norm ray packing](ABC_PRIMITIVE_NORM_RAY_PACKING.md)
· [Quadratic reference-ideal split](ABC_MINIMAL_REPAIR_IDEAL_CONTENT_SPLIT.md)

The unrestricted repaired proof obligation can be stated
as an average modular incidence of support records with
bounded-height representatives. The true occurrence count
has the same exponent as the primitive quotient count. At
prime size \(Q>N^{2+\eta}\), the principal and bounded-order
multiplicative modes now give a power saving. The entire
remaining obligation is a saving for the signed Jacobi-sum
correlation with all three high-order transforms.
Explicit primitive dictionaries on every fixed balanced
cap profile show that absolute values can attain the
full elementary upper bound even at those large primes.
Summing one character variable exactly before taking
absolute values still leaves a lower bound of order
\(N^{3/2-o(1)}\) on such examples.
[Absolute Jacobi sharpness](ABC_ABSOLUTE_JACOBI_SHARPNESS.md) Smaller
aggregate moduli \(Q=N^\beta\), \(1\le\beta<3/2\), attain
the generic false-positive exponent \(3-\beta\) on explicit
bounded-height dictionaries with all three actual radical caps.
[True support occurrence](ABC_MINIMAL_REPAIR_SUPPORT_OCCURRENCE.md)
· [Jacobi-sum reduction](ABC_MODULAR_JACOBI_INCIDENCE.md)
· [False incidences at smaller moduli](ABC_MODULAR_FALSE_INCIDENCES.md)

At the support primes themselves, minimally repaired
Pythagorean power families pass the first two layers with
\(Y^{1/M-o(1)}\) distinct target radicals. For polynomial
models the Wronskian bounds the total formal depth and
excludes simultaneous depth three. At actual primitive
integer parameters with fixed-degree forms, subpower
coefficients, and irreducible factors of degree at most nine,
passing depth three leaves at most \(H^{9/5+o(1)}\)
primitive parameter pairs. The extension uses a formally
low-depth factor and a uniform count for actual power
equations in its binary-form values. This leaves unrestricted integer occurrence open.
[Second-layer counterexample](ABC_MINIMAL_REPAIR_SUPPORT_OCCURRENCE.md#correct-heights-and-repair-through-the-second-support-layer)
· [Polynomial support depth](ABC_POLYNOMIAL_SUPPORT_DEPTH.md)
· [Primitive binary-form power counts](ABC_BINARY_POWER_TWISTS.md)

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


### Recent reductions for the balanced endpoint

At any fixed \(r<1/6\), rows with two distinct primes of
ambient occurrence at most \(X^{r\lambda}\) incur an absolute
exponent loss of order \(\lambda\), with the margin depending
on \(r\). A single rare prime on one letter forces the other
two letters to have nonatomic mass at least \(1/3-r\) each.
These conditions hold across the whole balanced endpoint.
The two-prime saving transfers to least-repaired primitive
quotients, with no count of extra common multiples. Typical
affordable divisor records also give many rows with the same
powered divisor, and with controlled additional support
intersection. Near one third a count for these actual
common-modulus fibres must still retain the complete radical
budget and the numerical height restrictions.
[Rare occurrences and common powerful moduli](ABC_RARE_PRIME_OCCURRENCE.md)

One quantitative target is to find, in any putative saturating
box, a prescribed root record \(\mathcal R\) and two letters
\(U,V\) with admissible conditional groups for which
\[
        2H(\mathcal R)+
          \sum_{G\subset U\text{ or }V} H(\boldsymbol x_G\mid\mathcal R)
                     \le (2/3-\eta)\lambda\log X
\]
for an absolute \(\eta>0\). The weighted dictionary bound
(GD4) would then contradict saturation. Small divisor
packets provide many candidate records, but the common-divisor
star bound alone supplies lower bounds on their multiplicity.
An upper bound for the actual three-letter incidence is needed
to force the stipulated conditional compression.
[Conditional dictionary inequality](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md#conditional-dictionary-inequality)
· [Common rows across a shift dictionary](ABC_AGGREGATE_SHIFT_DICTIONARY_STARS.md)

For fixed affordable fibres, the projective packing gap
forces a positive lower numerical weight on typical rational
output secants even after excluding a prescribed polynomial-height
prime pool. Consequently nonconstant rational charts of a degree
fixed before the height limit can produce height \(n^K\) outputs
on only \(n^{o(1)}\) rows in the small-budget exponent sense, for each
fixed \(K\). This closes the small-output branch of the pure-term
two-mixed-factor reduction, including charts with opposing
orientations and more than two slots. Counting enough original
tuples whose transferred outputs obey the requisite radical cap
is still open.
[Rational output height on affordable fibres](ABC_MONOMIAL_TRANSFER_HEIGHT_AND_FERMAT_GEOMETRY.md#a-uniform-small-output-bound-for-rational-charts)

The Darmon–Granville argument uses full triangle-cover indices
and a fixed coefficient support. In contrast, the ordinary
weighted projective curve has the effective common orders as
its inertia, including genus zero for pairwise coprime high
orders. The exact common multiplier then lifts every primitive
ratio to an exact-power point, with its height depending on
the primitive valuations. For least repair, a count must
control the cost of these lifts at varying supports. The
fine valuation grid supplies such control on its sector,
but the small width of a valuation window alone imposes no
common exact divisibility.
[Effective orders and weighted Fermat curves](ABC_EXACT_POWER_REPAIR_COEFFICIENTS.md#the-effective-orders-determine-the-weighted-fermat-curve)
· [Fine valuation divisibility](ABC_FINE_VALUATION_DIVISIBILITY.md)
· [Darmon–Granville, Theorem 2](https://www.math.mcgill.ca/darmon/pub/Articles/Research/12.Granville/pub12.pdf)

In the #943 lifting calculation, let two edges have a common
right kernel \(d\) and normalized left values \(A_iS\), with
\(S=\alpha^kg^j\). Write their right values as
\(\kappa_2d^jc_i^k\) and allow targets \(N_i\) to vary.
Eliminating \(S\) gives
\[
 \kappa_2d^j(A_2c_1^k-A_1c_2^k)=A_2N_1-A_1N_2 .
\]
The lifting step transfers after fixing the nonzero weighted
target difference, since the right-hand side is precisely
the Pell or Thue datum. Fixing an ordinary shift \(N_1-N_2\)
does not fix that datum when \(A_1\ne A_2\). A useful aggregate
transfer would need control of the weighted differences across
the actual target support dictionary, followed by recovery of
the common scale from a target or from an equally strong
arithmetic datum.
[Pell recovery](ERDOS_943_PELL_RECOVERED_KERNEL_GCD_WEIGHT.md#2-the-actual-cubic-and-norm-equations)
· [Mixed-power lifting](ERDOS_943_MIXED_POWER_KERNEL_GCD_WEIGHTS.md#2-an-exact-fixed-level-lifting-problem)


### Fixed shifts, high fullness, and fresh radicals in secants

The fixed-value strategy can dispense with a radical budget
on the prescribed shift. Above half the cheaper dictionary
cost, fixed-residue moments force both variable summands
into bounded lists of narrow roots with wider cofactors.
Each original root of positive numerical weight costs at
least the normalized fibre exponent. Fixing both complements
and releasing one high root per side gives one third of
the combined dictionary cost in original height shells.
The nondegenerate monomial estimate releases the whole
bounded list and gives
\[
       D\le\max\{\min(W_U,W_V)/2,B\}+o(\lambda)
\]
for bounded high groups, where \(B\) is the cost of the
low-degree singletons. The fixed coefficients can have
positive numerical height for this last bound.
[Two-root fixed-shift count](ABC_FIXED_SHIFT_TWO_ROOT_BOUND.md)
· [Arbitrary shifts with high free degrees](ABC_FULL_SHIFT_DICTIONARIES_AND_CODEGREES.md#arbitrary-shifts-with-high-free-degrees)

For each \(\eta>0\) and sufficiently large fixed \(k\),
the resulting ordered count of sums of two \(k\)-full
integers at every prescribed \(N\) is
\[
                      O_{k,\eta,\varepsilon}
                         (N^{(1/2+\eta)/k+\varepsilon}).
\]
The proof allows arbitrary common divisors. It fixes the
parts of the primitive quotients supported on the common
divisor and preserves all other original high valuations.
The unequal-fullness version uses half the smaller
reciprocal fullness when the two thresholds are proportional
to \(k\). The same fixed-shift count controls codegrees
after eliminating the shared row factor in a varying-shift
graph, using a full shift divisor or the high-degree
hypothesis on both targets.
[Pointwise powerful sums](ABC_FIXED_SHIFT_TWO_ROOT_BOUND.md#sums-of-two-highly-full-integers-at-every-target)
· [Weighted codegrees](ABC_FULL_SHIFT_DICTIONARIES_AND_CODEGREES.md#weighted-differences-in-a-varying-shift-graph)

At saturation, a typical joint divisor record of product \(G\)
inside an allowed group has occurrence frequency
\(G^{-1}X^{o(\lambda)}\). This gives many rows with extra
overlap arbitrarily small relative to their common record
inside the group. The bound across unrestricted positions
is anisotropic, and an explicit abstract support model
attains the remaining loss. In a maximal balanced fixed-value
fibre, the two free letters each need at least \(1/6\)
nonatomic mass, with the atom at scaled degree zero
included in the atomic budget.
[Joint divisor occurrences](ABC_DIVISOR_RECORD_RECIPROCITY.md)

A near-full family of \(m\)-tuples from an affordable fibre
of exponent \(r\) requires primitive radical cap at least
\(4r\) under a polynomial transfer. For a signed product
of input ratios, the two pure products recover the original
rows after discarding costly overlaps. If their charged
radical caps are \(A_1,A_2\), the remaining output term
needs fresh radical exponent at least
\(2mr-\min(A_1,A_2)\). Fresh primes avoid all three original
letters of each input row. For balanced pairs in cheap
records the required fresh radical in the projective
secant approaches \(2r\). This excludes low fresh-radical
compression. The unrestricted count would need control
of many actual secants, including the sector with large
fresh radical.
[Polynomial thresholds and fresh radicals](ABC_POLYNOMIAL_TRANSFER_RADICAL_THRESHOLD.md)


### Dispersion of secant valuations and fixed-shift applications

Finite-window concentration of fresh secant powers is a
power-sparse relation in a positive-exponent affordable
fibre. Fix two anchors with a large common neighbourhood.
The ratio of their secant values is a nonconstant projective
function of the query. A cheap full divisor on one
original query letter fixes that ratio modulo the
divisor. When boundedly many shrinking windows capture
its two secants' numerical height, the wide-cofactor
bound recovers the projective query at arbitrarily small
remaining cost. This argument permits a bounded radical
defect on each edge. After projective peeling, it gives
an \(r/2\) lower bound, up to fixed allowances, on
exact-power defects in a fibre of exponent \(r\) in base
\(n=X^L\), for every order above a sufficiently large
fixed threshold, including orders chosen by individual edges.
The diffuse numerical degree law also has root cost at
least \(r/2\) when it has no mass at zero or infinity
and its root cost is integrable.
[Finite-window secants and defect counts](ABC_FINITE_WINDOW_SECANT_GRAPHS.md)

The half coefficient for highly full fixed-shift pairs has
a strict saving for each fixed pair of proportional fullness
thresholds. At equality the more highly full term would
exhaust its root-height budget at the minimum degree.
It can then be compressed into one narrow root and a
vanishing-height cofactor. Mutual modulus locks force
boundedly many original high positions on both terms,
and their simultaneous release contradicts the positive
fibre exponent. For arbitrary common divisors and fixed
coefficients, the underlying fixed-shift estimate gives
additive energy at most \(n^{5l/2+o(1)}\) for a complete-value
dictionary of exponent and raw cost \(l>0\) with no atom at
scaled degree zero. A single large high-degree packet also
reduces varying-shift codegrees, even if division by the
common gcd removes some of its primes.
[Strict fixed-target saving](ABC_FIXED_SHIFT_STRICT_HALF_SAVING.md)
· [Energy and large-packet codegrees](ABC_NONPRIMITIVE_FIXED_SHIFT_ENERGY.md)

The two-cheapest dictionary bound charges the quotient of
a progression gap by its fixed step. For endpoint radical
caps \(\rho_1,\rho_2\) in base \(X\) and parameter gap at
most \(X^\tau\), the exponent for all close pairs is at
most
\(\tfrac12\min(\rho_1+\rho_2,\rho_1+\tau,\rho_2+\tau)\)
up to the small-budget allowance. The step radical still
enters the hypothesis guaranteeing an anatomic budget.
For a fixed physical gap, that hypothesis can be dropped:
the large step remains a fixed coefficient on a short
letter. For highly full endpoints no shortness of the
physical gap is needed either. A bounded-degree release
can instead keep one determined low root on the shift.
An upper count of occupied blocks in unrestricted
affordable fibres remains open.
[Progression pair counts](ABC_POWERED_PROGRESSION_OCCUPANCY.md#discounting-the-step-in-a-two-letter-bound)
· [Arbitrary large coefficient and a determined low root](ABC_FIXED_SHIFT_STRICT_HALF_SAVING.md#a-prescribed-coefficient-of-arbitrary-height-with-one-low-root)


Even without a step radical, physical gap or fullness
condition, the occupancy of one unit progression block of
parameter width \(X^\tau\) with endpoint radical at most
\(X^\rho\) is bounded by the exponent
\(\min\{\tau,3\rho/5,(2\rho+\tau)/5\}\), up to the
small-budget and height allowances. The proof combines
the two-letter count with a low-degree error on only one
tall letter and a modular uniqueness count for the low
part when the step is large. It charges each prime that
becomes low after gcd division against the radical of
the other original endpoint.
[General-step progression occupancy](ABC_GENERAL_STEP_PROGRESSION_OCCUPANCY.md)

## Output profiles of secants and products

A near-full relation of primitive monomial ratios at the fresh radical
threshold produces a second normalized \(1/3\) family. The pure
output degree laws are explicit weighted dilations of the input
letter laws. For a cross secant, the output pair-excess inequalities
force fresh nonatomic root mass at least
\(2/3-\min_{U\ne V}(c_U+c_V)\) in the original \(X^L\) units.
At the minimum nonatomic total this is \(1/3\), and such secants
must have numerical height exponent 2 in the original height
variable, bounded scaled degrees and no escaped numerical weight.
The transfer applies at arbitrary
nonatomic totals, with the height and tail loss controlled by the
specific input pair excess.

Fixed binary polynomials of the products require at least \(m/3\)
fresh radical exponent per distinct irreducible factor, up to
arbitrarily small fixed losses outside a power-smaller set.
Nonlinear factors are counted by fixed-field norm recovery,
extending the Pell counting mechanism and the minimally repaired
cyclotomic split. An adaptive exact \(q\)-th power quotient has the
stronger floor \((m/3-o(1))\tau(q)\) for its secant. A macroscopic
degree grid cannot carry the nonatomic input mass, so any power
order on a near-full relation has \(qL\to0\). This permits using
the small-budget counting bounds after taking the exact roots.
[Primitive secant transfer and adaptive power-class bounds](ABC_SECANT_ENDPOINT_TRANSFER.md)

For two fixed columns of the arithmetic graph, one compressed
valuation core suffices for the fixed-cofactor and reduced-ratio
codegree counts, with the other target's complete valuation list.
This gives a power deficit for multiplicatively consistent
rectangles in dense eligible graphs. A prime-labelled Latin
dictionary shows why marginal cofactor cardinalities alone cannot
force a ratio deficit. The unresolved step is to derive a
cofactor-ratio deficit or a comparable constraint on small-radical
secants from the original ternary equation at arbitrary profiles.
[One-sided compressed-core codegrees](ABC_COMPRESSED_CORE_WEIGHTED_CODEGREES.md)

## Radical growth under unequal powers and repair defects

For a family of primitive triples with two tracked letters and a
third tracked up to a cofactor of vanishing numerical height, the
normalized counting exponent is strictly less than the combined
radical mass at positive scales on the first two letters. This
prebalance inequality excludes original balanced boxes with too
much weight at low orders on two letters. It also rules out
bounded recorded radical in products of balanced input ratios
when at least half of their individual weights are sufficiently
small relative to the largest. In the product result the fresh
radical floor can exceed any fixed multiple of the input
budget. The associated primitive target must have a positive
numerical contribution from prime powers below any prescribed
positive scaled order, with no prior cap on its total radical.
[Prebalance count and imbalanced product transfer](ABC_IMBALANCED_PRODUCT_RADICAL_GROWTH.md)

A prescribed power class of a ratio of two independent affine
targets gives a fixed nonzero binary norm equation. Uniform
recovery in that norm bounds adaptive off-grid radical defects
on graph edges and rectangles. In a dense affine graph the
least target repairs also have large rectangle defects when the
fullness order is sufficiently high. This excludes multiplicative
separation of the target repair classes even if the separating
row and column factors have arbitrary height. The primitive
quotient formulation removes arbitrary common-multiplier copies.
[Affine power defects and least target repairs](ABC_RECTANGLE_POWER_CLASS_DEFECT.md)

For several norm forms of independent input directions, the
resultants from each omitted row and the discrete Loomis--Whitney
inequality fix the radical shared between the forms at negligible
exponent cost. Each fixed norm value then recovers the entire
tuple up to a subpower factor. With \(m\) slots and \(s\) forms
this forces joint fresh radical exponent at least \(sm/3\),
in the limiting convention and under the stated resultant
hypotheses.
[Joint multivariable norm bound](ABC_TENSOR_NORM_FRESH_RADICALS.md)

## Moving affine counts and arithmetic graph configurations

The fixed-shift count now gives coefficient \(2/9\) on any
moving noncoordinate affine line of actual primitive triples,
with all three radical caps retained. The prebalance inequality
handles the case in which support stripping leaves numerically
small quotients. When the weighted shift is relatively short,
the coefficient is \(1/6\) if the cheapest letter is available.
This gives a uniform \(2/9\) bound for primitive moving-offset
families \(m,m+\kappa,2m+\kappa\), and a \(1/6\) bound for
fixed \(\kappa\), up to arbitrarily small coefficient allowances.
It strengthens the small-budget unit-offset result of Zeng and
covers his fixed nonunit-offset case.
[Affine slice counts and equality profiles](ABC_AFFINE_SLICE_RADICAL_COUNTS.md)

Near-quadratic linear tripartite graphs can avoid balanced
six-cycles by Kong, Shangguan and Wang's construction. In
contrast, exact addition on complete labels forces
\(\gg E^{10}/N^{15}\) eight-edge trades for \(E\ge K N^{7/4}\)
and part cap \(N\). Their three sets of letter labels are
simultaneous arithmetic parallelograms. This supplies a
fixed-size configuration with a count at the endpoint
exponent using actual addition.
[Balanced cancellations in actual additive graphs](ABC_BALANCED_LATIN_TRADE_OBSTRUCTION.md)

With no affine label conversion, relative odds on common
columns give adaptive defect counts on three targets for
orders at least three, and on four targets for every order
at least two. A near-full graph with sufficiently high-order
least target repairs separating modulo powers is therefore
excluded directly in the outside projection.
[Homogeneous power stars and least target repairs](ABC_ADAPTIVE_HOMOGENEOUS_POWER_STARS.md)

For external modular tests, an additive four-shift estimate
for Fermat curve projections and two arithmetic differencing
steps extend the distinct-prime false-incidence obstruction to
total test exponent \(3/2\le\beta<38/25\). These boxes have no
integer equality and the three representatives are fixed across
the tested primes. The result constrains a proposed modular
count below the exponent needed by the true-occurrence test.
[Three-prime modular false incidences](ABC_THREE_PRIME_FERMAT_NOISE.md)

## Fixed targets, joint valuation classes and simultaneous families

A fixed complete letter gives half the smaller of the other
two radical caps as a counting exponent, also after subtracting
radical supported on any prescribed polynomial-height prime
pool. At balanced caps the complete-value fibre exponent is
at most \(1/6\), throughout the degree profiles. A complete
letter dictionary below exponent \(1/6\) therefore saves from
\(1/3\). For fixed-letter fibres the exponent is strictly
below both free letters' positive-scale masses. A free letter's
escaped high-order weight or a short-height gap also gives a
saving from the opposite free letter's nonatomic bound.
[Fixed-letter estimates and spectral restrictions](ABC_FIXED_LETTER_RADICAL_CAPS.md)

For least target repair of order \(M\), a single target radical
\(R\) supports at most \(R^{1/2+o(1)}\) primitive quotients at
balanced caps. The exact exponent allowance depends on the
balance tolerance. Thus a near-full family at repaired height
\(Y\) has target radical dictionary exponent at least
\((1/2-o(1))/M\). Several weighted valuation-class constraints
also combine through the least common multiple of their effective
orders. After removing one reference support and any necessary
high repair tails, a coefficient-uniform generalized Fermat
bound handles all three resulting powers. This gives a subpower
normalized fibre count even when individual effective orders
are small and the resulting fine grids have diffuse limits.
[Target radical fibres and simultaneous power classes](ABC_FIXED_LETTER_RADICAL_CAPS.md#size-of-a-minimally-repaired-target-radical-fibre)
· [Prior fine-grid sector](ABC_FINE_VALUATION_DIVISIBILITY.md)

Near-quadratic additive complete-label graphs contain simultaneous
binary shift allocations between any fixed list of pairs of
letters. Each allocation is an original triple. On a fixed fibre
the distinct shifted values have pairwise disjoint supports outside
a resultant pool, and the fixed-letter estimate shows that this
pool has negligible intersection with each radical at the balanced
complete-label exponent. A Behrend construction gives primitive
additive graphs at the same density that exclude the three
placements of a common shift.
[Translation cubes and fresh radicals](ABC_ADDITIVE_TRANSLATION_CUBES.md)

Pastén's arithmetic derivatives give a lattice inequality across
an arbitrary family of primitive rows. Simultaneous divisibility
can be added at each prime and valuation depth through the rank of
the two other letters' valuation differences. The common fixed
moduli alone have a capacity bound. Additional rank correlations
across rows or valuation depths are needed to turn the inequality
into a count.
[Arithmetic derivative family inequalities](ABC_ARITHMETIC_DERIVATIVE_FAMILY_LATTICES.md)

## Popular target supports and low-complexity complete letters

A modular characterization of popular target radicals uses
independent prime-test blocks. To test multiplicity at scale
\(N^u\), each block can have total prime exponent just above
\(2-u\), up to the balance tolerance, with a fixed number of
repetitions. Thus the upper half-budget fibre thresholds can
be tested at block exponent just above \(3/2\). At the critical
exponent the bounded-height Fermat boxes have polynomially many
false popular targets. The fixed-target Fourier completion gives
an explicit obstruction in the three-prime range as well as
an obstruction to using several undersized blocks.
[Modular tests for target multiplicity](ABC_MODULAR_TARGET_POPULARITY.md)

Pointwise bounds on fixed values and prescribed prime pools do
not imply the target multiplicity estimate on their own. A
finite-field evaluation construction gives an abstract graph
with total exponent one in repaired base \(N\), target count
exponent \(1-u\), and multiplicity exponent \(u\), for every
limiting \(0\le u\le1/2\). Its encoded integer labels obey the
balanced height and radical bounds and the pointwise pool
discount. They can pass any prescribed fixed modulus while
having at most one actual integer sum on their edges.
[Pointwise target multiplicity obstruction](ABC_TARGET_MULTIPLICITY_POINTWISE_OBSTRUCTION.md)

A bounded-rank affine model with sufficiently small maximal
integer minors on any one complete letter gives a uniform
saving at all balanced profiles. This includes progression
boxes of small positive numerical volume and approximate
polynomial descriptions with fixed degree and dimension.
The proof uses rank descent after cheap full-power locks.
When the designated letter is wholly atomic, both opposite
letters are diffuse and affordable divisors first lift their
complete values to the same affine model.
[All-profile bounds for small feature models](ABC_COMPLETE_LETTER_SMALL_FEATURES.md)

One complete letter or pair ratio in a multiplicative coset
of any fixed polylogarithmic rank also gives an explicit
saving at all profiles. A valuation basis separates the
values, and its prime partition sum allows a discount on
the original three radical caps. Separating bases can give
a count even without a rank cap. Three disjoint separating
prime sets for any one letter or pair ratio give coefficient
\(2/7\), by combining their entropy bounds with the other
two radical caps. Simultaneous packings on all three
letters give a stronger bound.
[Valuation cosets and packing counts](ABC_VALUATION_SPAN_COSETS.md)

## Short real windows and approximate valuation recovery

At every prescribed small coefficient \(\eta>0\), a sufficiently
short power interval at the numerical height contains at most
\(X^{\eta L+\varepsilon}\) integers of radical at most \(X^{AL}\),
for fixed \(A\) and all sufficiently small fixed \(L\). The proof
either records a cheap powered divisor larger than the interval
or counts the bounded list of large primes in a real strip.
This gives the corresponding local estimate for high-full
integers, and it makes coarse bins of each tall complete letter
lose arbitrarily little normalized entropy.
[Short intervals with small radicals](ABC_SHORT_INTERVAL_RADICAL_FILTER.md)

Combining this local estimate with modular rejection shows that,
after imposing the short real near-equality window, each
independent prime block for target multiplicity \(N^u\) can have
total exponent just above \(1-u\) in the balanced normalization.
At that same scale one can discard nonunit witnesses with
power-small loss in the true popularity tail, and remove every
fixed list of exceptional multiplicative character modes
pointwise. Complete-residue signed second moments are obstructed
on unfiltered prime-product dictionaries. A sufficiently small
signed high-character contribution on the actual filtered target
dictionary remains to be proved.
[Real-filtered popularity and signed character modes](ABC_MODULAR_POPULARITY_NONPRINCIPAL_MODES.md#unit-tests-at-the-smaller-real-filtered-prime-scale)

Valuation projections may now be used with explicit conditional
entropy losses. A projection missing only logarithmic valuation
rank incurs vanishing normalized entropy cost. In particular
three disjoint projections with such small kernels retain the
\(2/7\) coefficient. An abstract polynomial evaluation model
with balanced integer labels satisfies occurrence-frequency
constraints and still obstructs extracting inexpensive
projections from those data alone.
[Incomplete valuation projections](ABC_APPROXIMATE_VALUATION_SEPARATORS.md)

Long additive growth gives a further criterion. A bounded
normalized \(k\)-fold sumset for one complete dictionary at
\(k=\lfloor n^u\rfloor\), \(u>0\) fixed, gives a global
saving via Nguyen–Vu's inverse theorem. Short fixed-rank
progressions modulo sufficiently many smaller primes give
another saving via the Croot–Mao–Yip gluing theorem. An exact
primitive digit graph with subpower doubling shows why
fixed-order energies and primitivity alone cannot supply the
small-minor feature descriptions. It has the endpoint height
and cardinality scales and imposes no small radical cap.
[Long sumsets, modular gluing and the digit obstruction](ABC_HIGH_DIMENSIONAL_ADDITIVE_FEATURE_OBSTRUCTION.md)

## Two-base near-equality fibres and quantitative test obstructions

On tall shells, weighted sums or differences of two powers with
short bases have \(X^{\eta L+\varepsilon}\) near-equality fibres for
any fixed \(\eta>0\) and sufficiently small fixed \(L\), when at least
one order exceeds a fixed threshold. The coefficients and target may
move with height, and the orders may grow polylogarithmically.
The contact calculation uses a Wronskian in a single normalized
arithmetic parameter. Its degree is bounded independently of the
orders. For near-proportional opposite terms, either a rational
cross determinant fixes the reduced power ratio or the monomial
Wronskian is Vandermonde. This controls the real-window false pairs
on actual power-record fibres in any omitted-letter orientation.
If the available complete coefficient records per target number
\(N^{f+o(1)}\) in base \(N=X^L\), the sufficient block exponent is
\(\beta>f-u\) with a fixed margin at multiplicity \(N^u\).
The strict saving for short primitive summands already excludes
vanishing numerical height at saturation. The missing input is a
power discount on the coefficient records in arbitrary diffuse
valuation profiles.
[Weighted power near-equalities](ABC_WEIGHTED_POWER_NEAR_EQUALITY.md)
· [Strict saving for a short summand](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-strict-saving-when-the-short-term-is-small)

Reduced fractions with small joint radical are sparse in projective
windows of relative width \(X^{-(h_1+h_2)+\theta}\) at numerator
and denominator heights \(X^{h_1}, X^{h_2}\). Recording a cheap
full divisor of sufficiently large value makes each cross
determinant vanish in that window. The estimate can replace the
missing-rank hypothesis in an approximate valuation packing if
the projections localize the pair ratio to sufficiently few of
these numerical bins.
[Projective short-interval radical filter](ABC_PROJECTIVE_SHORT_INTERVAL_RADICAL_FILTER.md)

A bounded ramp on the real part of the high Jacobi modes, averaged
over independent primes, recovers the true target tail up to
constant threshold changes and a power-small error at the
filtered scale \(\beta>1-u\) with fixed slack. The actual third
target dictionary and its least-repair height are retained.
A uniform polynomial approximation to the ramp on the whole
multiplicity range requires power-growing degree for
\(u<1/2\). A full interpolation model attains all intervening
target multiplicity exponents while using every summand support
once and respecting the pointwise pool discounts. Its integer
encodings need not satisfy the real window.
[Capped signed target popularity](ABC_CAPPED_SIGNED_POPULARITY.md)

On comparable ratios, every continuum real filter formed by
absolutely summing separated logarithmic characters costs
\(\gg X^{(1-\theta)/2}\) in coefficient variation for window
\(X^\theta\) at height \(X\). The order is sharp. This rules
out a transfer that takes only absolute sums of separate
frequency bounds with fixed savings in base \(N=X^L\).
In another direction, the formal radical of a binary-form
model bounds its real contact with the sum equation. For
pairwise coprime fixed-degree formally high-full forms, the
real window admits only boundedly many primitive parameter
pairs per shell, even with subpower moving coefficients.
These obstructions leave the actual nonprincipal correlations
on arbitrary radical-constrained dictionaries open.
[Mellin cost of the real filter](ABC_MELLIN_REAL_FILTER_OBSTRUCTION.md)
· [Polynomial real-filter obstruction](ABC_POLYNOMIAL_REAL_FILTER_OBSTRUCTIONS.md)

## One-letter power grids, rich shifts and target moments

One high valuation grid bounds repeated differences of its complete
letter. In an actual addition graph, Cauchy–Schwarz on pairs through
a shared vertex retains the other two complete dictionaries. With
equal caps \(N^v\) on those dictionaries and defect radical at most
\(N^f\) on the grid, the resulting exponent is at most \(3v/2+f\)
in base \(N=X^L\), up to arbitrary allowances. At the minimum
nonatomic sector this gives \(1/4+f\). The cyclic anchor reduction
also gives a quantitative neighbourhood of that sector. A version
for a short real window uses the signed weighted-power fibre count
for long differences and short-interval radical sparsity for short
ones.
[One power grid and additive differences](ABC_ONE_POWER_GRID_ADDITIVE_GRAPHS.md)

In a near-full complete-label graph at the minimum nonatomic scale,
rich common-target shifts capped by \(\operatorname{rad}(h)\le N^b\)
have exponent at most \(\min\{b,1/18+b/3\}\), up to arbitrarily
small allowances. The proof transfers each shift to many tall
primitive difference triples and pays for the discarded common
multipliers by their supported radical. In particular a full
\(N^{1/6-o(1)}\) shift dictionary needs radical at least
\(N^{1/3-o(1)}\). At that critical radical cap, the derived
triples saturate a pair inequality with two copies of the
same original valuation profile. Their shifts must be full-height,
with bounded scaled support and no escaped numerical weight.
This leaves the possibility of rich shifts with radicals beyond
the critical cap.
[Radical amplification for rich shifts](ABC_POPULAR_SHIFT_RADICAL_AMPLIFICATION.md)

On minimally repaired primitive quotients with opposite even
exact powers of fixed order \(q\), a derived sum of the two
linear cyclotomic factors and one original root improves the
bound across the central repair band. The coefficient relative
to \(1/M\) is at most \(2/(\tau(q)+1)\) in the limit of small
budget and balance tolerances. When \(3\mid q\) or \(4\mid q\),
another derived primitive triple uses only radicals from the
cyclotomic factors of the intervening difference. This gives
coefficient \(1/\tau(q)\) relative to \(1/M\) at balanced caps.
Both bounds use the unique least multiplier of each quotient.
[Cyclotomic factor identities and minimal repair](ABC_MINIMAL_REPAIR_CYCLOTOMIC_SPLIT.md#additive-identities-for-opposite-even-powers)

In the pointwise model, every fixed finite list of moments of
target multiplicity is compatible with power-separated tails at
any prescribed exponent below \(1/2\). The two graphs can even
share their integer dictionaries and the same exact moment values.
Both obey the pointwise prime-pool constraints and keep total
edge count at the critical exponent. This specifies the ambiguity
left by a finite list of moments when used to test an intermediate
tail without the real relation.
[Matched moments in the pointwise graph model](ABC_CAPPED_SIGNED_POPULARITY.md#identical-finite-moments-with-different-critical-tails)

## Quadratic norm fibres and polynomial descriptions of rich shifts

A full divisor of one square root gives a compact ray class constraint
on the exact valuation layers of the third letter. The resulting
estimate works in both real and imaginary quadratic fields without a
numerical height-balance assumption. It excludes balanced saturation
whenever any pair has combined odd-valuation radical mass below
\(2/15\), with a fixed margin. The cost of the remaining class torsion
can also be weighted by which primes divide the actual layer indices.
[Weighted square ray class moments](ABC_WEIGHTED_SQUARE_RAY_CLASS_MOMENTS.md)

For square summands with an exact power target, Gaussian factorization
gives additional bounds before minimal repair. At balanced caps the
original-height coefficients are \(1/6\) for a cube target and \(1/9\)
for a fourth-power target, up to allowances. For weighted squares
and an odd exact target order, the same count uses the part of the
order recovered by principalizing the norm ideal. With full ideal class
variation, the alternative projective bound combines the already paid
summand kernels, target defect and any additional full-prime records.
At cubic target order the average \(3\)-torsion bound gives a strict
saving in the diffuse summand sector \(2f+3g<1/3\), where \(f\)
is the summand parity-kernel budget and \(g\) is the cubic-defect
radical budget. These bounds count original primitive quotients with
their unique minimum repair height.
[Gaussian power target splitting](ABC_GAUSSIAN_POWER_TARGET_SPLITTING.md)
· [Weighted norm descent and projective recovery](ABC_WEIGHTED_GAUSSIAN_CYCLOTOMIC_DESCENT.md)

For rich shifts, a polynomial description \(h=t_s+sP_s(z)\) with
\(\deg P_s\ge2\) and a cheap coefficient dictionary is excluded
independently of the radical of \(h\). A full-size rich-shift dictionary
requires coefficient exponent at least \(1/6\). The same threshold
holds near high powers when the coefficient divides the exact shift.
There are primitive dense graphs with tall prime rich shifts after
dropping the complete-letter radical caps. The large-radical shift
problem thus retains the radical constraints on the original letters.
[Polynomial descriptions of rich shifts](ABC_RICH_SHIFT_POWER_CARRIERS.md)

Jain's smooth-summand theorem gives a decomposition of each
sufficiently large integer with a polynomial largest-prime
cap on one summand and a polylogarithmic largest-prime cap
on the other. The hypotheses bound the largest primes.
They do not bound the product of distinct primes of the three
complete letters at a fixed small polynomial budget.
[Jain, Theorem 1](https://arxiv.org/html/2611.21136v1#S1)

A common multiplier twists the Frey curve by a scalar quadratic
character. The projective torsion at every level is unchanged,
including its inertia at primes away from the level and \(2\).
At level \(2q\) its ramification support there is the primitive
\(q\)-power defect. Making the repaired target a high power cannot
supply missing primitive level-lowering hypotheses.
[Frey inertia under target repair](ABC_FREY_REPAIR_PROJECTIVE_INERTIA.md)

## Mixed power-root ratios and higher-degree norm fibres

Passing to the ratio of two power roots splits the branch section
for the remaining letter into \(d\) sections. For fixed power
kernels, this gives height-free bounds exponential in the number
of defect primes for mixed signatures such as \((6,2,9)\),
\((6,6,3)\), and \((3,3,5)\). The estimate also allows moving kernels
with finitely many prime-to-\(p\) ratio classes. With fully moving
kernels the count is reduced to a weighted class-group moment.
For the height-free \((6,6,3)\) estimate, the variable class groups
can be taken to be two quadratic \(3\)-torsion groups. The Euclidean
\((3,3,3)\) estimate already gives the bounded-height defect-budget
bound for fully moving kernels in this signature. For each fixed
\((6,2,3^a)\), \(a\ge2\), a cyclic isogeny reduces the moving-kernel
cost to a higher moment of the same quadratic torsion. This moment
would extend the bound even at finite height, where the Euclidean
estimate does not apply. Support caps below the one-third threshold
give corresponding bounds on the primitive quotients with their
least repair.
[Root-ratio mixed covers](ABC_ROOT_RATIO_MIXED_COVERS.md)
· [Euclidean grid point counts](ABC_EUCLIDEAN_GRID_POINT_COUNTS.md#primitive-quotients-on-two-euclidean-signatures)

A bounded-degree compositum preserves each base torsion subgroup
up to a bounded factor. Even requiring common completely split
rational primes cannot itself compress odd-order class records.
Those primes can realize independently prescribed classes in the
moving quadratic fields.
[Composita and class-record obstruction](ABC_ROOT_RATIO_MIXED_COVERS.md#bounded-degree-composita-and-class-records)

Higher-degree projective ray class moments use a full divisor on
a root, or on the binary form value when its relevant irreducible
factor has degree at least three. The latter case keeps the
unoccupied prime ideals over each modulus prime. Comparing
conjugates forces diagonal products and gives a direct
target-repair bound from divisors of the repaired target's
primitive value. At each fixed order \(q\), this excludes
balanced saturation when a pair has combined \(q\)-defect
radical mass below \(2/[3(q+3)]\), with a fixed margin.
An explicit Gaussian family rules out the corresponding
degree-two value-divisor bound.
[Projective ray class moments](ABC_PROJECTIVE_RAY_CLASS_NORM_MOMENTS.md)

Consecutive valuation degrees can pass every test in a fixed
finite list of local power equations. For the quadratic signature
there are actual primitive additive realizations with these
prescribed degrees and polynomial height. For most prescribed
supports, every bounded-height realization needs fresh radical
above \(X_0^{bL}\) for any fixed \(b<2\), where \(X_0^L\)
is the scale of the prescribed radical. This quantifies a
global constraint absent from the local tests. The exceptional
supports left by that argument can still carry a near-saturating
family.
[Simultaneous local power tests and fresh radical](ABC_SIMULTANEOUS_LOCAL_POWER_PROFILE_OBSTRUCTION.md)

## Joint free-letter compression and moving square kernels

On a balanced profile with a cheap full modulus from one letter,
Fourier moments retain joint root dictionaries on both other
letters. A bounded group, a group in a fixed number of shrinking
valuation windows, or one group of numerically short full prime
powers can carry a joint restriction. Every such partition needs
full marginal entropy at saturation, even above the minimum
nonatomic mass.

For fixed \(a\ge2\), the mixed signature \((6,2,3^a)\) has a strict
unconditional saving when the total signature-defect radical mass
is bounded below \(1/3\) by a fixed margin and either the full
odd-kernel product is numerically short or its radical is concentrated
at a bounded number of narrow degree windows. The count allows
fully moving kernels. A first-moment estimate for quadratic
\(3\)-torsion makes the high-torsion kernels a sparse joint dictionary.
In the residual sector of the mixed signature, saturation requires
a positive nonatomic mass on the odd valuation degrees of the
two square letters.
[Joint free-letter dictionaries and moving kernels](ABC_JOINT_FREE_LETTER_DICTIONARIES.md)

Near-full families and near-half fixed-letter fibres with their
stated individual caps contain a growing set of actual rows with
mutually private full prime powers and arbitrarily small radical
loss. The same selection works against every later polynomial-height
prime pool. The aggregate least common multiple can have
quasi-polynomial height, leaving open how to use all these powers
together in a count with the original numerical height restriction.
[Growing private prime packets](ABC_GROWING_PRIVATE_PRIME_PACKETS.md)

## Growing arithmetic tests with private prime support

At almost every edge of a near-minimum complete-label graph, one
can extract polylogarithmically many fixed-size rectangles rooted
at that edge. Every new complete value keeps a private packet and
all original radical caps. The new values are multiplicatively
independent modulo the anchor primes. Subpower thinning can still
destroy every prescribed growing biclique while preserving the
count exponent and all original caps.
[Rooted private rectangles](ABC_ROOTED_PRIVATE_RECTANGLES.md)

A growing list of affine differences in one actual fixed-letter
fibre carries growing fresh radical on typical queries. Outside
the listed determinants and anchor primes, distinct differences
have disjoint supports. Their weighted products carry power
defect proportional to a fractional matching of the effective
orders, simultaneously for all orders and weights. Valuations of
these private primes also force the degree of any field containing
all the corresponding power roots to grow, including the roots
of least-repair multipliers at sufficiently large repair order.
[Growing affine defects and radical extension degrees](ABC_GROWING_AFFINE_DEFECT_FIBRES.md)

The private packets supply projective tests and simultaneous
ray-class tests at all balanced profiles. A nonzero secant or
conjugate difference bounds the total weight of matching probes.
A sufficiently large fixed collection of joint ray tests makes
the ideal-product moment diagonal. Averaging the tests
individually still allows power-sized collisions in a primitive
Gaussian family.
[Projective and ray tests on private packets](ABC_PRIVATE_PACKET_RAY_SIEVE.md)

For the full odd-kernel residue moment, the uniform local
Davenport–Heilbronn-type bounds give rarity in sections
specified by small kernel divisors. An abstract dictionary
can have that section rarity and still lie in one residue
fibre. There is an additional unconditional regime when a
small omitted radical mass leaves a signed ratio of height
below the full modulus. The weighted residue moment on the
remaining diffuse odd roots requires further arithmetic.
[Conditional first moments and retained ratios](ABC_JOINT_FREE_LETTER_DICTIONARIES.md#conditional-first-moments-and-short-retained-ratios)

## Averaging full divisors and cross-form residue locks

For a fixed off-diagonal ideal-product pair, an arbitrary family
of overlapping full moduli has only subpower many passing
records. This gives an ambient/local moment bound which accepts
arbitrary record-dependent subdictionaries. Its savings on the
off-diagonal branches cannot by themselves exceed the record
cost. In particular the arithmetic mean alone cannot finish
the high-torsion odd-kernel estimate. A separate affine graph
count uses the three original radical caps and gives a strict
density loss when a fixed number of cheap full divisors have
combined numerical weight above the query interval.
[Overlapping modulus moments and affine codegrees](ABC_OVERLAPPING_MODULUS_MOMENTS.md)

A full divisor coming from a different coprime form can instead
provide the scalar ray condition on a fixed binary carrier.
At a balanced profile, consider weighted square summands with
joint parity-kernel mass at most \(f\), a cubic target with
defect mass at most \(g\), and one summand with parity kernel
mass at least \(k\). A strict saving follows for
\(5f/4+3g/2-k/2<1/6\), at sufficiently small fixed budgets.
There is no assumption on individual prime sizes. Choosing
the larger kernel gives the sector \(f+3g/2<1/6\).
For exact cube targets with highly asymmetric kernels the cap
approaches \(f=2/9\).
[Cross-form projective ray locks](ABC_CROSS_FORM_PROJECTIVE_RAY_LOCKS.md)

Private odd valuations give actual complete labels with
balanced quadratic signatures at almost every prime above
the square root of the individual radical scale. Burgess's
interval estimate controls the off-diagonal quadratic character
products. The selection works in a balanced near-saturating
family or a fixed-letter fibre under its near-half hypotheses.
The count controls a quadratic mode at relevant test primes
while the full signed three-letter test still requires joint
information.
[Private quadratic signatures](ABC_PRIVATE_QUADRATIC_SIGNATURES.md)

The affine defect argument applies directly to forms in both
free coordinates of a triple, by testing a pair on a role slice
where both forms move. The lower bound depends on that role's
complete dictionary cap. A query-dependent pool removes overlap
without losing any of the three original radical constraints.
In particular this controls primitive difference directions in
rooted rectangle completions and the degree of extensions
containing their least-repair roots. On the regularized
minimum-mass graph, a supplying opposite band above \(s>3\)
forces outside radical mass at least \(3/2-3/s\) on any
full-exponent class of primitive directions with bounded
raw radical and a fixed excluded prime pool.
[Two-variable affine defects](ABC_TWO_VARIABLE_AFFINE_DEFECTS.md)

For a three-point cover of arbitrary degree, the reduction
in normalized tame discriminant cost is bounded by the largest
prime-divisibility mass of the exact valuation degrees on each
role. On planted diffuse prime-degree profiles the discriminant
cost plus the excluded planted radical mass is at least
\(2/3-o(1)\) for genus at least one, uniformly in the cover.
These profiles have actual quadratic realizations. Their
uncontrolled fresh radical still prevents using them as
counterexamples to a global saving.
[Three-point discriminant costs](ABC_THREE_POINT_DISCRIMINANT_COSTS.md)

## Primitive affine secants and complete-letter dictionaries

For near-full balanced families, almost every pair of actual
triples has a primitive affine difference of positive numerical
height and radical at least \(X^{(1-\eta)L}\), for each fixed
\(\eta>0\) and all sufficiently small fixed \(L>0\).
Two numerically short directions to fixed anchors would give
bounded-rank homogeneous features of the entire original row.
A full prime-power divisor on one diffuse letter forces a
rank drop. Effective radical counts for tall directions and
the planar incidence theorem then give a pair inequality
at all balanced profiles. Any power-abundant class of
secants below the indicated radical threshold, with fraction
loss smaller than the gap in (SC3), forces a strict global
row-count saving. The same geometry with endpoint-conditioned
prime pools gives sequential fresh radical on every fixed
number of bounded-radical directions.

A separate normalized coordinate-pair count uses the two
pure-letter radical caps of its output. A physical
difference with radical at most \(X^{BL}\) for fixed \(B\)
on a pair relation of full exponent forces the corresponding
complete letter of both inputs to the minimum dictionary
exponent \(1/6\) in base \(X^L\). In the central
least-target-repair problem, a
critical multiplicity class with exponent \(u<1/2\) in base
\(Y^{1/M}\) must therefore have almost all distinct target
differences outside every such fixed radical cap. Bounded
primitive directions give the complete-dictionary conclusion
when the common scaling gcd also has bounded normalized
radical. The three original caps have yet to force any of
these bounded-difference regimes.
[Primitive secants and physical coordinate differences](ABC_PRIMITIVE_SECANT_RADICAL_AMPLIFICATION.md)

Fixing one endpoint of a physical shift also fixes a
complete letter of the normalized output. One can therefore
discount the prime support of any polynomial-height integer
chosen from that endpoint, with no raw cap on the shift.
In a complete fibre of exponent \(r\), each of any fixed
number of successive shifts must contribute at least
\(2r-\eta\) fresh radical mass, outside the first row and
earlier shifts, on all but a power-small fraction of tuples.
At intermediate target multiplicity this gives an adaptive
form of the distinct-target difference bound and a separate
bound for shifts within the same complete target. At maximum
multiplicity, full-exponent common-target pairs with bounded
raw shift radical force minimum dictionaries for both moving
roles as well as the target.
[Fixed-endpoint shift counts](ABC_FIXED_ENDPOINT_SHIFT_RADICAL_COUNTS.md)

Actual comparable primitive sums with only the summand and
target caps give wholly diffuse least-repaired target families
at full occurrence exponent whose target differences exceed
every fixed normalized radical cap on all but a power-small
fraction of pairs at late budgets. The cyclotomic norm
count identifies the third original radical cap as a required
arithmetic input for any general difference extraction.
[Two-cap obstruction](ABC_PRIMITIVE_SECANT_RADICAL_AMPLIFICATION.md#a-two-cap-obstruction-with-actual-diffuse-least-repairs)

## Moving curves, the scaling gcd and fixed-fibre triangles

A moving plane curve of any fixed degree contains at most
\(X^{(1/3-\delta)L+\varepsilon}\) balanced shell rows for
an absolute \(\delta>0\). The budget threshold can depend
on the degree. A cheap full divisor of a diffuse letter
reduces the curve equation to a one-variable congruence.
A single-letter residue moment gives the saving after
enumerating the lock. Interpolate a polynomial through
fixed query anchors and use the least common multiple of
their affine secant gcds to force exact vanishing.
This puts the common scaling gcd on an endpoint-dependent
polynomial-height pool up to a numerically short factor
for almost every ordered pair. The exponent saving in the
curve bound is independent of the degree fixed before the
height limit.
[Moving curve count and local gcd pools](ABC_BOUNDED_DEGREE_CURVES_AND_SECANT_GCD.md)

The two-letter dictionary inequality tolerates a third
output with a fixed coefficient, a bounded-radical recorded
factor and a sufficiently short unrecorded cofactor.
Combined with the local gcd pools, this removes the scaling
gcd hypothesis from the primitive secant coordinate test.
If a critical pair relation has bounded outside radical
on one primitive coordinate, the corresponding
second-endpoint complete dictionary is minimal.
Symmetric tests give both conclusions. Intermediate
least-repaired target multiplicity therefore forces
large outside radical on primitive target differences
even after peeling the gcd. The same count allows a short
unrecorded part of the primitive coordinate itself.
In particular it forces positive numerical weight at
low scaled orders outside the pool after removing all
sufficiently high-order factors. Independently, a critical
relation with numerically large gcds forces at least two
minimum complete dictionaries. A role whose escaped first
moment is below the numerical gcd threshold also becomes
minimal.
[Short-cofactor and primitive coordinate counts](ABC_UNRECORDED_COFACTOR_SECANT_TRANSFER.md)

In a fixed complete-letter fibre of any positive normalized
exponent \(r\), the two moving letters retain positive
nonatomic mass even after recording arbitrary
polynomial-height prime pools. This rules out short
normalized difference triangles on a power-large fraction
of triples of actual complete values. At a full-exponent
triangle relation with radical cap \(n^{6r+o(1)}\),
incidence counting yields an output family of exponent
\(2r-o(1)\). Its primitive additive triples attain the
balanced coefficient at the new radical scale, including
at intermediate least-repaired multiplicities.
[Fixed-fibre primitive triangles](ABC_FIXED_FIBRE_TRIANGLE_RADICALS.md)

## Exact scaling and mixed prime packing at minimum dictionaries

A full-exponent tagged relation of fixed-fibre triangles
in a minimum-dictionary family can be thinned to one exact
scaling gcd and three fixed original residue classes with
vanishing normalized loss. Rare physical gaps carry a
power-smaller set of tagged triangles. The surviving
normalized sides belong to minimum dictionaries globally.
Pinning one complete target and one triangle value then
gives an injective primitive-output relation. At the
critical cap the new endpoint has minimum dictionaries,
even if its total nonatomic mass is above the minimum.
Two of its complete labels map by the same fixed affine
map back into a single actual original fibre. In a single
original fibre of any positive exponent the exact scaling
selection still works. The critical pinned output reaches
maximum target multiplicity at the new radical scale.
[Exact scaling and pinned triangle transfer](ABC_HEREDITARY_MINIMUM_TRIANGLE_TRANSFER.md)

For all minimum-dictionary profiles, the compact upper
cutoff on scaled orders extends from the nonatomic
spectrum to actual prime powers. With all roles fixed in
residue classes modulo \(G\ge X^\gamma\), the finite
cutoff is \(\max\{6,12(1-\gamma)\}\) and the nonatomic
cutoff is \(6(1-\gamma)\). One can obtain a common
modulus of numerical exponent approaching the sum of
the three escaped moments by first fixing the high-order
tails and then aligning residues successively.
[Compact prime packing and escaped-moment alignment](ABC_MINIMUM_DICTIONARY_COMPACT_PRIME_PACKING.md)

The mixed star count uses both other roles for anchors.
It combines a repeated high-order prime on one letter
with arbitrarily small primes on another. Write
\(H=1-\gamma\) and let \(b\) be the largest upper support
point among the three nonatomic profiles. A prime of
limiting cost \(t\) and order \(s\) obeys
\(kst+b(1/6-kt)\le H\) for each integer \(k\) with
\(0<kt<1/6\). In particular any finite mass above
\(6H+\sqrt{36H^2-6Hb}\) must lie on a single atom of
cost \(1/6\) in its own role, exhausting that role's
atomic mass. This restricts the higher-order end of
the spectrum even when the nonatomic masses are
arbitrary at minimum complete dictionaries.
[Mixed packing across the original roles](ABC_CROSS_ROLE_COMPACT_PACKING.md)

Two obstructions identify the arithmetic still needed.
The carry-free primitive construction can attain minimum
triangle-ratio entropy on one large complete fibre while
avoiding bounded-rank small feature models. An explicit
two-cap construction with diffuse complete values can
also keep primitive sums and a prescribed repaired-style
target multiplicity from zero through the half-budget
scale, while almost every primitive target coordinate
has unbounded normalized radical. A fixed-target
construction has the same behavior for all three
normalized triangle sides. The omitted summand radical
exceeds every fixed budget power on almost all its rows.
[Minimum triangle-entropy model](ABC_CRITICAL_TRIANGLE_ENTROPY_MODEL.md)
· [Two-cap primitive difference obstruction](ABC_DIFFUSE_TWO_CAP_PRIMITIVE_SECANT_OBSTRUCTION.md)

## Numerical widths and additive relations among reciprocal coordinates

A complete-letter fibre of positive exponent retains a
uniform numerical width after fixing any affordable pair
of full divisors on its moving letters. Their product
is a single CRT modulus and the length to compare it with
is the smaller of the two moving heights. In particular
the sum of their escaped numerical weights has a uniform
gap below that height. The same gap is available after
buying a finite upper quantile, with a separate rounding
rule for the nonatomic mass. This applies to every fixed
positive complete-target multiplicity stratum, including
least-repaired primitive quotients. A complementary
three-letter test uses the actual prime at the rounding
threshold. If a fractional quantile beats the projective
height threshold, a surviving full-exponent row requires
one actual prime with simultaneously large raw cost and
numerical weight. Fragmenting a limiting atom among small
primes cannot supply this obstruction.
[Fixed-fibre widths and actual-prime quantiles](ABC_FIXED_FIBRE_VALUATION_WIDTH.md)

On any full critical triangle relation in a fixed fibre
of positive exponent \(r\), an opposite-vertex pin gives
the two minimum dictionaries required to test the remaining
side. Repeating this separately for each side bounds its
finite spectrum in terms of its height and the two other
escaped moments. On the original side, a one-letter
variable-mesh estimate makes a polynomially short
relative interval or a positive-height full unit modulus
cost a fixed share of the exponent. At minimum nonatomic
mass it forces full physical triangle heights and excludes
positive high-order weight in the exact scaling gcd.
[Two-dictionary packing and the original affine images](ABC_TWO_MINIMUM_DICTIONARIES_TAIL_PACKING.md)

On a fixed chain of critical minimum-dictionary transfers,
a terminal complete value has distinct affine lifts into
every prior level. If \(M_2,M_3\) record the full high tails
on the two unpinned parent legs, then
\(\operatorname{lcm}(M_2,M_3)\mid G_jG_{j+1}\) for the
two exact scaling gcds. The rank-one width bound caps the
product of all the scales uniformly in the fixed chain
length. These results constrain chains that remain
available on the retained outputs. The original three
radical caps still do not force critical triangle
relations or make them hereditary across every
refinement.
[Pinned chains and inherited full moduli](ABC_PINNED_TRIANGLE_CHAINS.md)

The reciprocal coordinate \(1/a+1/b-1/c\) has reduced
denominator \(abc\). A fixed exact linear relation among
distinct such coordinates requires every denominator
prime outside the coefficient and target pool to recur.
Anisotropic radical pricing charges each prime at its
first occurrence. It gives a power saving when a positive
share of the radical occurs at least three times, already
lies in the fixed pool, or leaves unequal first-occurrence
costs on the three roles. The argument also counts relations
of least-repaired coordinates, where the denominator
becomes \(Z_M(c)abc\) with unchanged prime support. The
unstratified additive moments have a generic rational
denominator proof. Their additional abc content is the
quantitative dependence on radical overlap and on the
three role costs.
[Reciprocal coordinates and radical overlap](ABC_THREE_CAP_RATIO_SUMSETS.md)
· [Bourgain–Dilworth–Ford–Konyagin–Kutzarova, Lemma 5](https://denka.web.illinois.edu/rip_stoc.pdf)
· [Bloom–Maynard, §1](https://arxiv.org/pdf/2011.13266)

For the projective coordinate \(F(z)=z-1/z+1/(1-z)\),
a full power \(p^e\) on the last row using that prime
supplies a homogeneous congruence modulo \(p^{2e}\).
Fixing the high tails at small radical cost lets these
congruences determine a primitive row at an intermediate
step. At every fixed odd order there is a uniform exponent
loss for exact linear relations with that many distinct
active rows. At even orders the macroscopic exclusively
repeated prime packets of surviving tuples form a perfect
matching. Nontrivial complete-letter Latin trades admit
a separate \(B/9\) exponent saving at normalized total
radical cap \(B\), including in the repaired dictionaries. A complete primitive
additive graph constructed by CRT has precisely those
reciprocal relations at any fixed order. The construction
has no small radical budget at the required counting scale.
The new estimate leaves the problem of deriving a
quantitative reciprocal compression, abundant odd relations
among distinct rows, or nonmatching recurrence from the
joint radical and height constraints.
[Reciprocal compression and the primitive CRT example](ABC_THREE_CAP_RATIO_SUMSETS.md)

For this route, a sufficient additional estimate is the following
fixed-order statement. Find an odd \(k\ge3\) and a fixed coefficient
height bound such that every family within some fixed normalized
distance of \(1/3\) yields a weighted relation in \(F(a/c)\)
with at least \(n^{k/6-\delta_k/2}\) solutions on distinct rows.
The constant \(\delta_k>0\) is the loss furnished by (RS8)
for that order and height bound. The weights and target may
be selected from the family after taking its height box,
and are fixed while counting its solutions. This would
contradict (RS8) at all sufficiently small fixed \(L\).
The known dictionary and fixed-fibre constraints do not give
this lower bound.
[Exact odd-relation gap](ABC_THREE_CAP_RATIO_SUMSETS.md#near-maximal-relations-pair-the-repeated-prime-supports)


The elementary quadratic norm recovery from Erdős #943
is uniform when both norm parameters are fixed. The
least-repaired quotient leaves the two parameters mobile
across the target supports. For example, writing
\(a=k_a x^2\), \(b=k_b y^2\) with squarefree kernels gives
\(c^2-k_a k_b(2xy)^2=(a-b)^2\), with varying norm and
discriminant. The small-fundamental-unit estimate for the
ordinary Pell equation gives a gain in the symmetric family
\((n,n+1,2n+1)\) with coefficient loss of order
\(\lambda^2\), through an exponent saving of order
\(\lambda^3\). The generalized Pell identity on an
arbitrary fibre does not put it on that ordinary
small-unit branch. The counting problem here needs a
bound for the joint occurrence of the moving records.
[Erdős #943, uniform elementary norm lemma](ERDOS_943_PELL_RECOVERED_KERNEL_GCD_WEIGHT.md#1-a-uniform-elementary-norm-lemma)
· [Zeng, §§5–6](https://hexagonmath.org/pdf/2611.00611v1)
· [Diao, Theorem 1.3](https://arxiv.org/pdf/2610.28001v1)


## Exact pair recovery and moving reciprocal targets

The even pairing obstruction for exact scalar relations
in \(F(z)=z-1/z+1/(1-z)\) can be removed on sorted
primitive rows. A pair satisfying \(F(x)=rF(y)+s\)
maps to a genus-two cyclic cubic curve with four split
branch points over \(\mathbb Q(\sqrt{-3})\). The two
degenerate affine transformations force equality or
reflection. Outside them, split \(1-\zeta_3\) descent
bounds the rank in terms of the primes of a
polynomial-height branch discriminant. Uniform
Mordell–Lang then makes the last pair cost subpower
uniformly in the affine coefficients. Expose the pair
with the largest exclusive prime packet last. At any
fixed order \(k\ge3\), the normalized exponent for
distinct active rows improves from \(kB/6\) to
\(kB/6-B/[3(k-1)]\), with arbitrarily small allowances,
under total radical cap \(n^B\). This applies directly
to the least-repaired quotient coordinate, which
depends only on the primitive ratio. There is also a
scalar fixed-target version for arbitrary common gcd
in powerful representations.
[Genus-two recovery and all-order relations](ABC_PROJECTIVE_RECIPROCAL_GENUS_TWO_RECOVERY.md)
· [Gao–Ge–Kühne, Theorem 1.1′](https://arxiv.org/pdf/2105.15085v4)
· [Poonen–Schaefer, Proposition 6.2 and Lemma 12.9](https://math.mit.edu/~poonen/papers/descent.pdf)
· [Dimitrov–Gao–Habegger, Corollaries 1.3 and 6.1](https://arxiv.org/abs/1904.07268v2)

Joint reciprocal coordinates have an elementary
two-fraction pair recovery and the same sequential
saving. They give a sharper fibre count with prescribed
outside pools. For fixed powerful targets this
argument can be applied before restricting the common
gcd: prescribe its divisor class, use the free
radical caps outside the actual target, and expose
short primitive quotients separately. The rates also
hold for the projective three-coordinate version
unchanged by least repair. Component targets may
move in intervals when their reduced denominator sizes
permit bounded reconstruction from spacing.
[Joint reciprocal pair recovery and fixed powerful targets](ABC_JOINT_RECIPROCAL_PAIR_RECOVERY.md)
· [Browning–Elsholtz, Theorem 1 and §2](https://www.math.tugraz.at/~elsholtz/WWW/papers/papers33FINAL2013.pdf)

Small target denominator *radical* also has direct
moving-target bounds. For scalar \(F\) the singleton
outside prime mass raises the baseline exponent by
at most \(\mu/6\) at target radical cap \(n^\mu\).
The odd-order saving persists if this cost is small.
For two consecutive signed reciprocal powers one
obtains a loss at every fixed order, globally and
on a fixed complete letter. At an exclusively shared
prime, the two powers force the signed complete
values to agree modulo a higher power. At a last
prime they recover the signed letter from the ratio
of the two known prefix sums. Fixing a positive
compact valuation band on the tracked roles lets the
argument use either alternative at an intermediate
row. Six separate role sums extend this result to the
actual least-repaired letters, even when the least
multiplier varies. There is a complementary exact
support estimate using a longer block of consecutive
powers. A local Vandermonde annihilator bounds the
full divisor of each active label missed by all the
output denominators. Under the compact-band condition,
small denominator radical in that longer block forces
exact label cancellation.
[Moving scalar output](ABC_THREE_CAP_RATIO_SUMSETS.md#odd-relations-with-variable-output-denominator)
· [Consecutive powers and actual minimal repair](ABC_RECIPROCAL_POWERS_VARIABLE_TARGET.md)
· [Vandermonde support inequality](ABC_RECIPROCAL_VANDERMONDE_SUPPORT.md)
· [Bloom–Kuperberg, Theorem 3](https://arxiv.org/abs/2312.09021)

The cyclic determinants of a projective pair give
a second primitive additive triple, with the full
powers at exclusive cancellation primes surviving
squared. Eisenstein factorization recovers the
ordered input pair from the exact determinants at
subpower cost. A small *radical* of the three
cofactors, together with two sufficiently full cyclic
packets, allows the two-full-letter theorem to exclude
that sector at a fixed exponent loss. The recovery
pays for the possible common factor of the
determinants, and short primitive outputs are counted
separately. The remaining input here is simultaneous
control of the numerator radical of the projective
difference and of the cyclic allocation of high
orders.
[Cyclic determinant transfer to two full letters](ABC_CYCLIC_PROJECTIVE_PAIR_TRANSFER.md)
· [Two-full-letter saving from endpoint diffusion](POWERFUL_TRIPLES_ABC_ENDPOINT_TRANSFER.md#corollary-a-small-budget-saving-with-two-full-letters)

For the exact scalar relation route, one can now ask
for near-maximal abundance at any fixed order
\(k\ge3\) on distinct rows in a near-one-third family.
For example, at normalized budget \(B=1\), a fixed
polynomial-height weighted relation with
\(n^{k/6-1/[6(k-1)]}\) distinct-row solutions would
contradict the new pair bound. For the moving-target
route the corresponding input is enough tuples with
small outside radical in the joint reciprocal
denominators. The graph density and cardinality
results in hand supply neither lower bound, so the
uniform saving for the unrestricted triple count
remains open.


## Radical growth of projective sums and forced trade defects

For a fixed tuple of \(k\) capped rows, the private primes of
the individual denominators survive in a shifted projective
sum. Recording the output numerator and exposing a numerically
larger block of the rows lets those private full powers
determine the product of the unseen denominators by CRT.
Consequently a numerator cap below
\(\lfloor(k-1)/2\rfloor/3\) in base \(n\) gives a fixed
power loss from the full tuple exponent \(k/3\). There is
also a joint bound for two shifted numerators with fixed
rational weights and shifts. They and their difference
give a primitive additive output, while its denominator
recovers the input tuple at small cost. The two-cheapest-letter
inequality forces fresh numerator radical mass at least
\(2k/3-o(1)\) on the pair, in addition to denominator
radical mass \(k-o(1)\). Taking shifts \(0,-1\) gives
total output radical mass at least \(5k/3-o(1)\).
The argument covers the projective coordinate \(F\) and
the Frey \(j\)-invariant. Its fixed-fibre version uses
the half-cap price on the moving letters. For an actual
fixed powerful target it allows arbitrary gcd and the
cost of excluding the target support from both shifted
numerators is independent of the order.
[Projective numerator recovery and output radicals](ABC_PROJECTIVE_SUM_NUMERATOR_RADICALS.md)

One-row fixed-degree norm forms give an associated
fresh-radical threshold. An affordable full divisor from
a positive nonatomic original letter places the norm
value in a projective ray class. A norm field with
discriminant cap \(X^{gL}\) then requires normalized
fresh radical at least \((2/3-g/2)_+\) at saturation.
Pairwise resultants separate distinct forms away from a
fixed affordable prime pool. In particular the Eisenstein
norm \(a^2+ab+b^2\) and the cyclic cubic numerator of
\(F(a/c)\) together require fresh mass \(4/3-o(1)\).
[Half price for fixed-degree norm radicals](ABC_PROJECTIVE_SUM_NUMERATOR_RADICALS.md#half-price-for-fresh-radicals-of-fixed-degree-norm-forms)

The supersaturated eight-edge trades in actual
minimum-dictionary graphs do not give scalar reciprocal
identities. Their signed defect has an explicit
nonzero sign. Under balanced radical caps the forced
trades give \(n^{5/6-o(1)}\) distinct denominator
radicals of mass at least \(4-o(1)\) in base \(n\).
The privacy of the twelve label radicals follows
by fixing four independent labels and discounting
overlap on one complete-letter fibre. Separate
dense primitive graphs at the endpoint height scale
have unique six-term \(F\)-sums and hence show the
limits of using small label sumsets alone.
[Forced trade defects and the endpoint-height example](ABC_PROJECTIVE_RECIPROCAL_GRAPH_OBSTRUCTION.md)

Two consecutive positive powers of the scalar projective
coordinate yield a saving linear in the number of distinct
active rows, at every fixed order. The exponent can be
chosen before the order because two-row cancellation
bounds the numerical weight carried by any one eligible
pair. In a random exposure order, a fixed fraction of the
rows have enough earlier full powers to reconstruct their
coordinate and still have a positive later radical charge
to remove. Both output denominators may move within a
small outside-radical cap. The coordinate is unchanged
by least repair, and the fixed-powerful-target version
allows arbitrary gcd.
[Consecutive scalar powers and random exposure](ABC_SCALAR_RECIPROCAL_POWER_AMPLIFICATION.md)
· [Bloom–Kuperberg, Theorem 3](https://arxiv.org/abs/2312.09021)

The split genus-two pair recovery also applies to
fractional-linear relations. A liftable branch-preserving
exception requires an exact cube determinant for a
primitive integer matrix. If such a map relates two
rows with large common radical, all of that common
support lies in the fixed pool of its bottom-left entry
and determinant. A small map dictionary then permits
effective radical pricing. Near-quadratic reciprocal
sumset and Möbius product expansion hold on the rational
image of \(F\) at polynomial height, complementing
field-level incidence bounds. A sufficiently small
dictionary of fractional-linear maps for the
unlabelled balanced families has not been obtained.
[Fractional-linear recovery and arithmetic expansion](ABC_FRACTIONAL_LINEAR_PROJECTIVE_RECOVERY.md)
· [Zhi, Theorem 1.1](https://arxiv.org/html/2611.26661v1#S1)

These results make the low-denominator incidence
obligation more restrictive. The forced short trades
carry large output radical, while low-radical actual
sums are sparse in the full family. A general
improvement requires a mechanism producing some
different arithmetic configuration with retained count
exponent, or an upper bound that uses the proven fresh
radical growth to directly constrain the original
exceptional-set cardinality. The currently available
anatomic and generalized Fermat estimates supply the
one-row and fixed-fibre inputs for these transfers.
[Bernert–Browning–Lichtman–Teräväinen, Proposition 2.1](https://arxiv.org/html/2410.12234v2)
· [Browning–Verzobio, Theorem 1.2](https://arxiv.org/html/2608.24512v1#S1)



## Adaptive valuation grids and centered norm constraints

The fine-grid count extends to adaptive choices of \(m\) exact
orders across the three primitive letters. If every selected
order is at least \(K\) and the joint radical of primes served
by none of their own letter’s orders is at most \(X^f\), the
exponent is at most \(f+C_m K^{-m^2/(m^2-1)}+\varepsilon\).
The number of orders can increase sufficiently slowly as the
budget decreases. The proof freezes very large powered values
and releases the remaining bounded array of disjoint roots
in their actual ternary relation. The choices of \(m\) and
\(K\) precede the height limit. A prime-index valuation law
gives the same diffuse weak limit as a fine exact grid. It
evades every sublinear dictionary of exact divisibility orders.
Brun–Titchmarsh and rank-two compression show that it also
evades signed grids of bounded rank and coefficient-width
product \(o(1/\lambda)\). The law can be realized on coprime
integer letters. Its realization in an additive endpoint
family is still an open condition.
[Adaptive grids and the prime-index obstruction](ABC_FINE_VALUATION_DIVISIBILITY.md#a-bounded-number-of-adaptive-valuation-grids)
· [Prime degree indices in sparse signed grids](ABC_FINE_VALUATION_DIVISIBILITY.md#prime-degree-indices-also-evade-sparse-signed-grids)

For \(Q=a^2+ab+b^2\), any exponent-saturating subfamily of
the original count has \(\operatorname{rad}(Q)\) of normalized
mass at least \(2/3-o(1)\). More precisely, if the original
cap is \(c^\lambda\) and \(\operatorname{rad}(Q)\le c^{B_Q\lambda}\),
the limiting count coefficient is at most
\[
 \min\!\left(\frac13,B_Q,\max\!\left(\frac14,\frac{B_Q}2\right)\right).
\]
The proof splits according to whether any original linear
form supplies a cheap full divisor of positive numerical
height. The Eisenstein norm receives a half price in the
projective ray class in the first case. In the other case
finitely many exact valuation coordinates carry the full
heights, all three forms are tall, and bounded-root release
gives the \(1/4\) coefficient. A joint cap on the radicals
of \(abc\), its centered reciprocal companion, and \(Q\)
gives coefficient \(1/8\). The companion also gives a
six-branch cyclic twist. Its defect for every fixed prime
order \(\ell\) yields a height-free support-only count
\(O_{\ell,\varepsilon}(Y^{1+\varepsilon})\) when the joint
\(\ell\)-defect radical is at most \(Y\).
[Centered reciprocal companion and norm counts](ABC_CENTERED_RECIPROCAL_NORM_COUNTS.md)

## Polynomial relations on actual least-repaired heights

Write \(H(c)=W(c)\operatorname{rad}(c)^M\) with
\(M=\lceil K/L\rceil\). The support of the valuation
tail \(W(c)\) costs at most \(1/K\) in base \(n=X^L\).
In a polynomial relation on \(k\) actual repaired heights,
split the radical roots by exact prime membership across
the \(k\) slots. A packet whose membership vector has
the same dot product with every monomial exponent vector
only dilates the height tuple. Every other packet is
pairwise coprime to the others and has exponent spread
at least \(M\), so coefficient-uniform nondegenerate power
release gives \(n^{k/K+\eta}\) records up to a height
subpower. If the invisible membership vectors span
rank \(r\), then \(r\) height coordinates fix a tuple
within each record. On primitive rows the corresponding
count is \(T^r F^{k-r} n^{k/K+\eta}\), where \(T\) is the
family size and \(F\) is its largest repaired-height fibre,
with the height-subpower allowance understood.
This yields a fixed loss from the full tuple exponent
across the balanced endpoint for every fixed relation
with \(r<k\), provided \(K\) is large enough.

For homogeneous linear relations the invisible span has
rank one. If \(D\) is the number of distinct repaired heights,
off-diagonal equal \(h\)-sums have count at most
\(D^{h-1}n^{2h/K+\eta}\) up to a height-subpower allowance.
The resulting \(h\)-fold sumset has asymptotically the
largest size allowed by permutation symmetry when
\(D\ge n^u\) and \(K\) is sufficiently large in terms
of the fixed \(h,u\). With a nonzero inhomogeneous
output there are no invisible packets. In particular,
for fixed weights and pool, the number of pairs with
nonzero weighted difference of outside radical at
most \(n^b\) is at most \(n^{b+2/K+\eta}\), again
with the subpower allowance. These statements use
the unique least repair. Primitive multiplicity is
separate and is controlled by the fixed-letter bound.
[Polynomial relations on least-repaired heights](ABC_REPAIRED_HEIGHT_POLYNOMIAL_RELATIONS.md)

A two-stage Pell and Thue transfer gives a different use
of actual repaired differences. A nonzero weighted target
difference, or one nonzero nonlinear factor of it,
recovers the two repaired targets. Recording the reduced
left ratio and the right power-class defect then recovers
the original pair of primitive rows. The common factor
of the right summands divides the actual recovered
difference, so no unrestricted common multiplier enters.
A joint outside radical cap \(A\) on the record yields
\(O(A Y^\varepsilon)\) pairs with repaired heights at most
\(Y\). In the small-budget formulation a record of
normalized mass \(s\), selected from a pool dictionary
of exponent \(\rho\), gives repaired exponent at most
\((s+\rho+2/K)/(M(1/3-\xi))+\varepsilon\).
It is strictly below \(1/M\) if \(s+\rho+2/K<1/3-\xi\).
Generic disjoint pairs can already pay too much on the
reduced left ratio, before counting the target datum.
Thus a use in the unrestricted count requires an
extraction of sufficiently many genuinely correlated
primitive rows with inexpensive coefficients and
difference records.
[Least-repaired weighted target differences](ABC_REPAIRED_WEIGHTED_TARGET_DIFFERENCES.md)
· [Mixed-power kernel-GCD lifting](ERDOS_943_MIXED_POWER_KERNEL_GCD_WEIGHTS.md#2-an-exact-fixed-level-lifting-problem)

A three-role arithmetic construction shows why it matters to
retain each unrecorded cofactor’s radical. It gives a primitive
matching at exponent \(1/3\) with diffuse exact valuation slots
on all three roles. Every union of slots attains the expected
record entropy, including affordable records of numerical
weight greater than one. The full marked numerical weight can
approach two. One letter has its actual balanced radical cap.
Both other letters have large unrecorded radicals, carried
entirely by bounded valuations, so they disappear in every
compact positive scaled-degree test. The residual radical of
each of the three complete letters is essential to the sought
incidence inequality.
[Three-role exact divisor record obstruction](ABC_THREE_ROLE_DIVISOR_RECORD_OBSTRUCTION.md)

## Independent repair and the additive graph

Separately repair each original letter and retain the fact
that the three resulting prime supports in a row are disjoint.
A prime membership pattern may therefore select at most
one coordinate of each row. For a polynomial relation on
observables built from the repaired values, only these
patterns can contribute to its invisible valuation space.
An observable that uses two distinct repaired letters in
one row already fixes both original supports and recovers
the row at subpower cost. This replaces the fixed-target
multiplicity charge in the corresponding relation bound.
Every prescribed bounded-degree projective curve meets
the independent-repair image in at most \(n^{3/K+\eta}\)
primitive rows up to a height-subpower factor. Collinear
triples of distinct repaired points cost at most
\(Tn^{9/K+\eta}\) in a family of \(T\) original rows.
A moving rational sum can also be counted using only its
numerator when every negative local exponent minimum is
unique. After recording small shared support, this gives
fresh numerator radical growth for sums of reciprocal
repaired products \(1/H_M(abc)\) throughout the balanced
endpoint.
[Prime membership and projective incidence](ABC_REPAIRED_PRIMITIVE_PROJECTIVE_INCIDENCE.md)

For degree \(d\) interpolation, take a minimally dependent
tuple of \(j\le\binom{d+2}{2}\) repaired evaluation vectors
with distinct two-coordinate ratios. Keep all maximal
minors and split each into its minimal vanishing subsums.
A private prime packet invisible to all these subsums
could vary one coordinate while preserving dependence.
That would force every degree-\(d\) form through the
other points to contain a particular coordinate line.
Replacing that factor produces a form that still passes
through those points and does not contain the line.
This contradiction accounts for every private packet.
After paying \(\rho\) for shared radical mass, the resulting
tuple count is \(n^{\rho+3j/K+\eta}\) up to height subpowers.
The same argument releases packets meeting several rows
when the remaining level of the linear system has enough
dimensions. A conic of six actual separately repaired
primitive rows exhibits an invisible mixed-column packet
outside this range. Screening larger overlap gives a
fixed saving on all fixed-order repaired interpolation
dependencies at the original one-third tuple scale.
[Simultaneous power release for repaired interpolation](ABC_REPAIRED_VERONESE_DEPENDENCIES.md)

The independent-repair modulus method also treats transforms
\(\mu_M(T)^j T^k\) with bounded \(j\ge1,k\ge0\).
An affordable prime packet of one original letter becomes
a very tall modulus after repair. The ratio of the other
two transforms is an exact high power times a wider
cofactor, whose numerical height stays bounded in the
original scale. This proves a strict gap across all
profiles if a small dictionary of coefficients gives
an actual ternary relation of these transforms. Canonical
inverse-repair weights do not supply the dictionary:
even modulo simultaneous scaling and coordinatewise
rational \(M\)-th powers, they recover the original
prime supports outside the high tails and hence have
at least \(n^{1/3-3/K-o(1)}\) classes at saturation.
At the minimum complete-dictionary scale, the forced
eight-edge trade typically has a nonzero repair defect.
Its canonical correction across four targets also has
near-maximal coefficient-class count.
[Moduli and trade defects for independent repairs](ABC_INDEPENDENT_REPAIR_RELATIONS.md)

For any fixed positive primitive additive pattern, a
universal sparse identity in the repairs and the original
labels must split into groups carrying the same repair
monomial, when the order \(M\) is large relative to the
number of terms and coefficient degrees. To see this,
construct a primitive affine pencil on which the letters
are squarefree outside a fixed pool. Each repair becomes
a constant times a high power of a distinct linear form.
The function-field unit bound then distinguishes the
high-power profiles in a minimal vanishing subsum.
The same pencil prevents a small coefficient dictionary
from extracting such a nondegenerate relation on a
full-exponent subset using the additive equations alone.

The need to retain the third radical also holds for
actual graphs with two genuine diffuse radical caps.
A dense Cartesian projection on those two letters has
two repaired-height dictionaries with nearly diagonal
fixed-order energies. It admits a matching of exponent
\(1/3\) with all three complete-label dictionaries at
the critical size. The third letter has outside radical
mass approaching \(2/3\) in base \(n\). This realizes
the ceiling of cardinality-only entropy and product
convolution bounds on primitive rows.
[Sparse universal repair identities](ABC_SPARSE_REPAIR_GRAPH_IDENTITIES.md)
· [Two complete caps and the third radical](ABC_TWO_CAP_REPAIRED_ADDITIVE_OBSTRUCTION.md)

## Path congruences and original divisor fibres

Target divisor records overlap across a labelled additive graph.
At a given prime and depth, the congruences identify the signed
summands in each connected component. The exact projected modulus
for a new vertex is the least common multiple of its path moduli
to already specified vertices. The corresponding lattice index is
the product of the prime powers supplied by maximum spanning forests,
one per prime. This counts fixed overlapping records at the cost of
their union radical and allows several paths to reconstruct a vertex
even when no common-neighbor comparison can do so.
[Path projections and affordable record counts (GL1–GL4)](ABC_GRAPH_PATH_DIVISOR_LOCKS.md)

The same components yield a global inequality on an actual
three-letter graph. Put \(n=X^\lambda\), select prime-depth
pairs \(\mathcal J\) with \(p\le n^\beta\) outside the fixed divisor
primes, and write \(\mathcal W_{\mathcal J}(t)\) for the sum of
\(\log p\) over selected \(p^j\mid t\). Write \(D_{\mathcal J}\)
for the sum of \(\mathcal W_{\mathcal J}\) over distinct
unordered cofactor pairs on each of the two summand sides.
For \(h>2\) and targets of polynomial height,
\[
 (1-2/h)\sum_{(i,j)\in E}\mathcal W_{\mathcal J}(C_{ij})
 \le D_{\mathcal J}
       +O\!\left(h(N_A+N_B)n^\beta\log X\right).
\]
Components whose two sides differ in size by a factor of
at least \(n^\delta\ge2\) contribute at most
\(4n^{-\delta}D_{\mathcal J}\) to the target sum.
In minimum complete dictionaries, where
\(N_A,N_B=n^{1/6+o(1)}\) and \(|E|=n^{1/3-o(1)}\),
the error and imbalanced component contributions are
power smaller than \(|E|\log X\) for \(\beta<1/6\)
and \(h=n^\theta\) with \(0<\theta<1/6-\beta\).
The general error is power smaller wherever the edge
exponent exceeds the larger summand dictionary exponent
by more than \(\beta+\theta\).

In any such sector with a lower bound
\(\sum_E\mathcal W_{\mathcal J}(C_{ij})\ge
\eta|E|\log X\) for fixed \(\eta>0\), a uniform estimate
\(D_{\mathcal J}=O(n^{-\sigma}|E|\log X)\) for some
\(\sigma>0\) would give a contradiction. Here the
prime-depth tests and the numerical lower bound must
be established for the actual third letter. For a
uniform count one must supply a fixed margin in the
edge and dictionary exponents, and impose the estimate
for all sufficiently small fixed \(\lambda\), with
\(X\) tending to infinity first. Large balanced components
can still occupy a small fraction of both dictionaries,
and the upper bound on their summed difference moment is
open. This criterion applies to the specified degree
and dictionary ranges.
[Prime-depth difference inequality (GL5–GL7)](ABC_GRAPH_PATH_DIVISOR_LOCKS.md#a-difference-moment-inequality-for-the-actual-edge-layers)

Inside a fixed two-divisor fibre, the target projection
preserves the entropy and collision probability of sampled
sums to exponentially small error up to a fixed power
below the square of the numerical lock margin. The square
scale is sharp for primitive rows with unrestricted radicals.
A separate projective packing argument recovers a tuple
from a normalized target template if the residual common
target content fits below the lock margin. Together with
full-support peeling, this gives a count priced by the
relative shape and the target-only valuation defect, at
arbitrary relative indices and with no cost for summand
defects. The needed inexpensive shapes and target defects
still require an extraction theorem.
[Weighted projection and target-only templates (PE1–PE2)](ABC_AFFORDABLE_FIBRE_PROJECTION_ENTROPY.md)

A polynomial-height pool capturing the numerical parts of
the two summands also controls a large common divisor of two
targets. The elementary count against a fixed divisor
\(M\mid c\), at outside-pool numerical height at most
\(\tau\log X\) on \(a/b\), is
\(O_\varepsilon(X^\varepsilon(1+X^\tau/M))\).
The anchored form permits a different prescribed pool at
each anchor. At sufficiently low outside height, the
quantitative Subspace Theorem additionally covers those
pairs by \(X^{o(1)}\) plane curves of uniformly bounded
degree in the original ratios. A Cartesian grid gives an
explicit counterexample to the polynomial dependence on
logarithmic coefficient height in the stated moving-coefficient
bounds of Noubissie v1. The fixed-shift result uses Evertse's
theorem directly.
[Almost-unit common numerator gcd bounds (AP1–AP5)](ABC_ALMOST_UNIT_COMMON_NUMERATOR_GCD.md)

With full divisors on all three roles, an extra projective
congruence of cyclic index \(J\) gives a split cubic with
explicit pair determinants and discriminant \((FGHJ^3)^2\).
The three residual letters are exactly the three linear
factors of a split order's index form. Coordinate primitivity,
the exact span of the original rows, and the conductor of
the split cubic order are explicit. An example at actual
height realizes a number of congruence classes of order \(J\)
while \(\operatorname{rad}(J)\) stays fixed and the first
letter satisfies its cap. It has large unlocked radicals
on the other two letters. The remaining occupancy problem
keeps the separate caps on all three factors of the index
form.
[Split cubic locks and index forms (SL1–SL6)](ABC_SPLIT_CUBIC_LOCK_LATTICES.md)


## Deep residue blocks at marked compact primes

In a near-minimum complete-dictionary graph, mark a
target prime \(p\in[n^\beta,n^{\beta+\gamma}]\), where
\(n=X^L\), \(0<\beta\le1/6\) is fixed, and a power \(P_p\)
of that prime dividing the target has numerical height
at least \(X^w\) for fixed \(w>0\).
Give each prime its empirical price from the distinct
complete target values. The effective radical count
uses all three original caps to rule out target
occupancy exceeding \(n^{1/6-\beta}\) by a fixed
power. In the two summand dictionaries, an oversized
residue class would contribute disproportionately to
the pair-collision budget, since two distinct values
can agree at only \(O_w(1)\) of these deep moduli.
Small or sparse blocks also lose a power by direct
edge counts.

With arbitrarily small fixed exponent tolerances,
almost all marked edges of a near-full graph have
three actual block dictionaries of size
\(n^{1/6-\beta+O(\eta)}\) and block edge count at
least \(n^{1/3-2\beta-3\eta}\). There are
\(n^{\beta-O(\eta)}\) active primes. Each prime occurs on at least
\(n^{1/6-7\eta}\) values in each summand dictionary, in at most
\(n^{\beta+\kappa+\eta}\) unit classes modulo \(P_p\).
Also, at almost every retained edge both summand
endpoints are incident to marked constraints with
modulus product at least \(X^{n^{\beta-2\eta}}\).
This rules out a vertex-wise construction with
polynomial-size products of its incident moduli.
The size of the product alone does not bound the
least representative of an actual residue vector.
[Deep residue inverse and endpoint products (RB1–RB3)](ABC_DEEP_RESIDUE_BLOCK_INVERSE.md)

A fixed number of such moduli can now be chosen
simultaneously on a large subset of either summand
dictionary. If that summand has a representation
\(A=\rho t^d\), with \(t\) of rational height at most
\(X^u\) and \(\rho\) from a dictionary of exponent
\(\sigma\), this gives a strict saving whenever
\(kw>2u\) and \(\sigma+k\beta<1/6\) for a fixed
integer \(k\). The height condition is \(kw>u\)
for integer \(t\). One may choose \(d\) adaptively
from a finite list at each fixed budget. The same
conclusion holds for a prescribed list of nonconstant
rational maps after removing the primes of the fibre
discriminants. A direct full-modulus argument also treats bounded exact
multi-grids on one summand and arbitrary dictionary
sizes in the opposite role. The count exponent
there is at most \(b+\sigma+\beta+\gamma\), up to
arbitrarily small allowances, when the opposite
dictionary has exponent \(b\) and the chosen
coefficients are units at their marks. No grid is
needed on the other two roles.
[Adaptive power representations and multi-grids (RB4–RB5)](ABC_DEEP_RESIDUE_BLOCK_INVERSE.md#a-uniform-saving-from-a-power-class-on-one-summand)

The three orientations can also share one numerical
difference budget when each prime is assigned to
one original role. This bounds the constant density
of a near-complete triangle graph with almost all
numerical weight in the tested layers. For fixed
mixed-role patterns, the local ranks at each prime
depth give a joint divisor of all difference minors.
A divisor exceeding the box bound forces an affine
rank drop. A separate application of the short-term abc theorem counts
summand values with short unit residue representatives
of controlled radical, when the modulus and its
cofactor have controlled radicals as well.
Block progressions give actual primitive examples with
a genuine target cap, diffuse tested depth packets,
and abundant difference moments. The same construction
admits generic large residue representatives and
violates the summand caps. Thus a general difference
moment bound needs arithmetic information beyond the
target cap.
[Three-role precision, short residues, and the target-cap example](ABC_CROSS_ROLE_PRECISION_AND_SHORT_RESIDUES.md)

An additive polynomial-evaluation model realizes the
marked block occupancies and bounded pair agreement
at numerical height \(X\), with primitive rows. Its
marks are formal. A direct interpolation argument
rules out realizing nearly all of them by actual
prime powers using a polynomial digit encoding of
degree small compared with the digit ranges. Exact
additivity forces affine labels on the full template
when targets are indexed by integer digit sums.
These obstructions use compatibility across layers
beyond separate pairwise numerical capacities.
[Marked incidence model and polynomial lift obstruction](ABC_DEEP_TARGET_LAYER_INCIDENCE_OBSTRUCTION.md)

In a split cubic lock lattice, the actual indices
on projectively related original rows have a joint
power-class count. If the product of the index
forms has enough branches after cancellation, one
fixed projective correspondence with outside defect
radical at most \(B\) gives \(O(BX^\varepsilon)\)
row tuples. For parity on two rows, this holds
with zero or one common branch. Two common branches
allow a genus-zero family of \(\gg X\) primitive
pairs in one parity class. The cubic power class
retains the necessary branches for every two-row
projective map. A fixed-target shift gives a
genus-one application, and a separate uniform Pell
count uses a fixed first row and the same joint
parity cap. These statements require an actual
correspondence and a bound on the joint defect.
[Joint index classes and fixed-target shifts](ABC_SPLIT_INDEX_POWER_CORRESPONDENCES.md)

The residue inverse turns the marked sector into a
counting problem for joint high lifts. A sufficient
estimate is a uniform power saving in the number of
band primes for which two small-radical dictionaries
with the simultaneous occupancies above can be matched
into the required capped target blocks. The block proof
uses positive fixed \(\beta,w\) and near-minimum complete
dictionary sizes. It yields no fixed exponent margin as
\(\beta\to0\), and it leaves the broad-profile endpoint
obligation open.

## Composite packets and projective counting

A composite packet of full original target powers can be
marked even when every prime has subpower size. Use one
packet per complete target and bin by the number of targets
assigned to the packet. If that number has exponent
\(1/6-s\), the packet dictionary has exponent \(s\).
The fixed-letter estimate gives \(s\ge\beta/2\) up to the
fixed tolerances, where \(\beta\) is the packet radical
exponent. The raw support count gives \(s\le\beta+\gamma\).
Same-side pair divisibility and the target cap then give
actual residue blocks at exponent \(1/6-s\) on each role.
A fixed list of \(k\) almost disjoint packets occurs at many
distinct original summand values. This yields a strict count
when \(ks<1/6\) and their combined numerical weight exceeds
one, with fixed margins for the small removed overlaps.
A short root or algebraic feature allows a smaller numerical
threshold. The residue reconstruction works through
finite algebraic correspondences, including nonrational
ones, by removing the primes at which the selected
feature is a multiple root.
[Composite packet inverse and simultaneous lists](ABC_AGGREGATE_DEEP_PACKET_INVERSE.md)
· [Algebraic features in the packet inverse](ABC_DEEP_PACKET_ALGEBRAIC_LIFTS.md)

The empirical packet price can be as large as the raw
radical price under one-letter hypotheses. A polynomial
evaluation code gives actual integer targets in dyadic
height shells with the required radical cap and diffuse
compact-order valuations on subpower primes. Fixing the
primes of a packet fixes enough code coordinates to
force price at least \(\beta\), even for adaptive marks.
The same dictionary avoids every prescribed polynomial-height
prime pool at fixed positive radical mass with a power
saving. A compulsory price discount therefore needs
information from the actual additive graph or its other
radical caps.
[One-letter obstruction to packet compression (AP5)](ABC_AGGREGATE_DEEP_PACKET_INVERSE.md#a-one-letter-obstruction-to-compulsory-packet-compression)

An absolute residue on a primitive input of a norm form
gives an ordinary Arakelov ray condition at its inactive
primes. A list of \(R\) allowed residues contributes
\(R^{1/2}\) to the count. This covers quadratic norm forms
with a deep target divisor, and variable fields of bounded
degree with bounded discriminant exponent in base \(n\).
Larger fixed moments absorb the class-number exponent
in the small-budget allowance. For a separable reducible
target form, one can fix any nonlinear irreducible factor.
At a prime where the packet divides another factor, every
prime above it in the chosen field is inactive. The modulus
therefore retains positive numerical weight even if the
packet divides linear factors only.

Together with the composite blocks, this bounds the
coefficient by \(1/6+f+(3s-\beta)/2\) up to fixed
allowances when the two summands are weighted \(q\)-th
powers with combined defect radical exponent \(f\).
The same ordinary norm count transfers to minimally
repaired primitive quotients with prescribed target-divisor
and absolute-residue dictionaries. It still requires
a bounded-degree norm representation or a related
affordable carrier.
[Absolute residues, composite norm packets, and primitive repair (AT1–AT7)](ABC_ABSOLUTE_TARGET_RAY_MOMENTS.md)


An independent obstruction concerns actual prime-power
realizations of the evaluation-layer model. Even with
arbitrary bounded-height integer labels and a fixed small
fraction of missing incidences, those prescribed layers
cannot have the stated moduli if every \(k\) of them have
least common multiple greater than \(2^kX\). Translations
along the level kernels make a \(k\)-fold mixed difference
divisible by the combined modulus, so it vanishes over the
integers. Stability of the difference identity on the cyclic
digit then forces a large constant class, contradicting
injectivity. The argument controls near-complete realizations
of this specific template. It does not rule out exponent-dense
subgraphs of vanishing relative density.
[Robust CRT obstruction for evaluation layers](ABC_EVALUATION_LAYER_CRT_OBSTRUCTION.md)

A separate count applies at arbitrary complete dictionary
sizes and throughout the balanced profile range. If an
image \(\Phi([a:c])\) of the original primitive ratio has a
projective description in boundedly many integer features
with small affine minors, that family saves below \(1/3\).
A full packet on one of \(a,b,c\) forces the image into its
corresponding branch class. This gives an exact affine
rank drop in the feature space, even at a pole and with
shared content in the two presenting coordinates. The
modulus comes from a diffuse band and need not contain any
individual prime of positive normalized size. The result
includes rational functions of boundedly many short
integer parameters, and it gives a strict count for
minimally repaired primitive quotients with such a feature
model on the original ratio. It does not extract a model
from the radical assumptions.
[Projective feature count and primitive repair (PF1–PF6)](ABC_PROJECTIVE_SMALL_FEATURES.md)

## Simultaneous residues and primitive quotient counts

Assign one full divisor to each original complete letter.
When the sum \(d\) of the three empirical prices is below
\(1/6\) by a fixed margin, simultaneous marking retains
near-full packet strata. Fixing the three packets \(F,G,M\)
and their three absolute residues leaves an exact affine
equation on intervals of width \(O(1+X/(FGM))\), with the
three separate outside radical caps. Each complete dictionary
in a surviving cell has exponent \(1/6-d\) up to tolerances.
For a fixed packet tuple, the occupied residue triples fill a
power-near-full part of the Cartesian product of the three
residue alphabets. They contain many complete boxes of any
fixed size when the three prices are positive. Each box gives
actual rows in every cell, with the overlaps specified at the
level of complete-letter residue classes.
A CRT construction realizes arbitrary fixed boxes with shared
primitive values, exact full divisors and the numerical
height when \(FGM\operatorname{rad}(FGM)\le X^{1-\theta}\)
for fixed \(\theta>0\). The three outside radical caps are
the missing ingredient in that construction.
[Simultaneous composite packet blocks (SB1–SB3)](ABC_SIMULTANEOUS_COMPOSITE_PACKET_BLOCKS.md)
· [Residue pair projections and boxes (SC1–SC3)](ABC_SIMULTANEOUS_RESIDUE_PRODUCTS.md)

At the minimum nonatomic mass, a purely diffuse packet has
empirical price equal to its radical exponent \(\beta\).
After pinning one opposite complete value and the packet,
the residual target in that original fixed fibre must retain
nonatomic mass at least \(1/6-s\). More generally the
price is at least \(D-\alpha-c_T+b\), where \(D\) is the
edge exponent, \(\alpha\) the opposite dictionary exponent,
\(c_T\) the original target's nonatomic mass, and \(b\)
the nonatomic mass removed by the packet. In the minimum
sector the numerical weight of a purely diffuse packet is
at most \(6s\). This precludes a direct height saving from
substituting lower empirical prices for radical prices there.
[Nonatomic packet price rigidity (RP1–RP2)](ABC_NONATOMIC_PACKET_PRICE_RIGIDITY.md)

For a generalized Fermat carrier, two absolute residue
lists can shorten both root progressions. A good full
part of a divisor of the third letter then improves the
surface cover before counting the resulting curves.
With orders \((16,16,2)\), index widths at most
\(n^{1/10+o(1)}\), a retained target divisor of fixed positive
numerical weight, and record cost at most \(1/4\), the
resulting small-budget coefficient is \(3/10\). Under the
stated lower target radical cap this transfers to a strict
count of minimally repaired primitive quotients. The proof
requires the two short root progressions in its input.
[Two shortened root progressions and a third divisor (DF1)](ABC_OPPOSITE_ROOT_LIST_DETERMINANT.md)

A separate fixed-target count treats every order \(q\ge2\).
Prescribing the \(q\)-power class of the other two letters'
product or ratio leaves only \(O_\varepsilon(X^\varepsilon)\)
primitive rows, uniformly in the complete target and order.
The proof allocates the class support to the original
coprime inputs, then counts a fixed value of a separable
binary form; at large orders it counts the short roots
directly. The conclusion also bounds pairs with a cheap
defect off an anchor support, including adaptively chosen
orders, and applies to odds unchanged by common repair.
The canonical quadratic norm lift also records twice the
target packet depth at inactive split primes. An explicit
family over distinct quadratic fields shows that a small
field dictionary cannot be inferred from the packet and
the outside norm radical alone. Both original summand
caps would have to enter such an inference.
[Canonical quadratic packet lift and uniform class separation (CQ1–CQ5)](ABC_CANONICAL_QUADRATIC_PACKET_FIELDS.md)

The determinant step uses Browning–Verzobio's affine surface
theorem with a side condition and their generalized Fermat
curve geometry. The additional count here runs the curve
projection in both shortened progression indices. The fixed
degree norm step in the fixed-target count is the standard
weighted binary-form argument. Its application here is
uniform in the power order and gives an outside-defect
pair bound in the actual packets. For comparison, the
recent \(0.584\) whole-value modulus certificate is a fixed
exponent uniform over \(0<\lambda\le1\); it supplies no
small-budget coefficient at this scale.
[Browning–Verzobio, affine surfaces with a side condition](https://arxiv.org/html/2408.11453v2#S1)
· [Browning–Verzobio, sums of three powerful numbers](https://arxiv.org/html/2608.24512v1#S1)
· [Zeng, fixed-degree norm argument and complete modulus certificate](https://hexagonmath.org/pdf/2611.00621v1)

## Shared-value arrays and arithmetic obstructions

A surviving simultaneous packet tuple forces actual difference
arrays of every fixed size \(m\). Each displayed letter depends
on two of the three indices. All \(m^3\) rows are original,
and a single complete value is shared across all rows with the
same pair of indices. This follows from grouping ordered
common-target tuples by their exact difference vectors, then
taking complete bipartite graphs in the grouped lists with
distinct residues. If \(d\) is the sum of the packet prices and
\(t=1/6-d\), there are \(n^{md+(3m-1)t-O_m(\eta+\kappa)}\)
ordered arrays. Conversely the \(m\) disjoint diagonal
triangles anchor a spanning tree. This bounds the number of
boxes by \(E^mK^{m-1}\), with \(E\) rows in the tuple and
\(K\) the maximum pair-class size. Thus a power saving for
arrays of one fixed size and a power saving for rows are
equivalent in the marked regime with positive price margins.
[Sharp shared-value extraction (LB1–LB4)](ABC_SHARED_VALUE_RESIDUE_BOXES.md)

For every prime order, most such boxes give an injective
valuation map on the full circulation space of dimension
\(3m^2-3m+1\). The stronger count treats all orders and weights
at once and forces a positive power of fresh radical in each
nonzero circulation product. It recovers a missing potential
from two nonproportional affine forms by a fixed-value binary
norm count. Enumerating a joint defect support once lets this
recovery propagate along a forest. More generally, to reduce
the rank of the valuation map to \(r\), the same support record
permits recovering any next potential with at least \(r+2\)
known neighbors. The resulting support cost increases with
the number of recoverable potentials.
[Cycle defect count (LB3)](ABC_SHARED_VALUE_RESIDUE_BOXES.md#cycle-power-classes-on-actual-boxes)
· [Joint defects and rank (BF1–BF2)](ABC_SHARED_BOX_FOREST_POWER_DEFECTS.md)

The one-letter graph count now applies to any order \(q\ge2\)
chosen separately at each value. Partition by the witness
order and use the fixed nonzero norm count on shifts within
each part. At minimum nonatomic mass, one original letter
with outside defect exponent \(f<1/12\) gives a strict
saving even when the captured orders are small. Within a
simultaneous packet tuple, the anisotropic shift count on
\(A\) saves if \(2f<t+\max(s_G,s_M)\) by a fixed margin.
The same conclusion holds in the other orientations. Small
fixed orders may still give a diffuse limiting profile.
[Arbitrary-order one-letter defect tests (OG6–OG8)](ABC_ONE_POWER_GRID_ADDITIVE_GRAPHS.md#an-arbitrary-order-defect-on-a-single-letter)

On a fixed complete-letter fibre, a small integer feature
model for a nonzero integral multiple of a free letter
forces a power bound through full-divisor rank descent.
This gives a secant inequality inside packet cells and
forces primitive direction radical greater than
\(n^{6t-\sigma}\) for almost all pairs. For two actual
shared-value arrays with the same residue lists, the \(m\)
disjoint diagonal cells jointly contribute fresh direction
radical greater than \(n^{m(6t-\sigma)}\) outside the entire
first array and a fixed normalized-height pool. In the same
generic pairs, every spanning tree in the first array
already pays the full \(md+(3m-1)t\) radical exponent up
to a margin. These support lower bounds add, and the local
discriminant inequality then charges both the tree and its
new primitive directions to any common powerful lift.
Both arrays retain their three original radical caps.
Using one seed triangle and then adding one potential at a
time gives a stronger bound. A known complete difference
fixes one primitive direction entry at divisor cost.
The seed contributes \(6t-\sigma\) and each of the \(3m-3\)
later directions contributes \(4t-\sigma\) in fresh radical
exponent outside the first array and previous directions.
The same fixed-fibre bound proves
uniform sparsity in short progressions of arbitrary step.
This excludes packet products sufficiently close to the
complete numerical height. The pooled valuation test also
bounds the sum of the three packets' numerical weights,
the escaped weight on two free roles, and their top residual
nonatomic \(t\)-moment by the smaller of those two heights.
[Local primitive secants, progression widths and pooled valuation tests](ABC_CELL_PRIMITIVE_SECANT_BOUNDS.md)

A ramified base change makes the three principal ideals
powerful without introducing a common ideal factor. Its
minimum root discriminant is determined, up to fixed-order
factors, by the required local ramification indices. For
one prescribed field the actual fractional ramification
weights yield the effective radical count, uniformly over
field degree. In the least-target-repair sector this bounds
the primitive quotients per field at an exponent tending
to \(1/(2M)\). For a shared-value box there is also a bound
over all variable fields at once. The root discriminant
and the prescribed order bounds control the unweighted
prime support of the spanning-tree labels. Comparing with
the forced boxes makes a common field of row-scale
discriminant power-exceptional when the box is large enough.
[Optimal ramification and field-cover bounds (RL1–RL8)](ABC_RAMIFIED_IDEAL_LIFT_COST.md)

## Private support and rigid transfers on shared-value boxes

At the shared-value box exponent, every original edge
between axes \(i,j\) has private radical exponent at
least \(2(t+\max(s_i,s_j))\) up to a fixed tolerance.
The counted primes are absent from all other original
edges. A fixed-anchor fibre also gives a private
nonatomic mass of at least \(t+\max(s_i,s_j)\).
The individual bounds can be added across all \(3m^2\)
edges. On synchronized pairs one can enforce privacy
across both input arrays. These estimates use a fixed
complete original letter and all original radical caps.
[Private prime support on actual boxes (PV1–PV5)](ABC_SHARED_BOX_PRIVATE_PRIME_SUPPORT.md)

Pairs of boxes can retain all residue matches or leave
some specified indices unmatched. A Latin family of
\(m^2\) actual cells partitions the \(3m^2\) cross edges.
For exactly matched residue lists, its primitive directions
jointly contribute fresh radical exceeding
\(n^{6m^2t-3m^2\sigma}\) outside both input arrays and a
fixed normalized-height pool. Each new packet survives
primitive reduction because it belongs to only one raw
difference edge. This improves the linear number of
directions in the seed-vertex count. If the numerical
weight \(w\) of the three common packets exceeds \(6d\)
by a fixed margin, almost no such direction stays within
the original relative budget. The tree-support costs
of both input arrays and the new directions also add
in a common-field root-discriminant bound. Within a
single array, fixed-label shifts give an independent
private radical floor on every capped difference
between potentials on the same axis. Their joint
outside support has a corresponding unconditional
lower bound.
[Private differences, Latin directions and same-axis shifts (PA1–PA4)](ABC_PRIVATE_RADICALS_ON_ARRAY_DIFFERENCES.md)

The full edge space, of dimension \(3m^2\), has independent
power classes on generic boxes for every order \(q\ge2\).
If two active weights meet, a single free potential can be
counted through a fixed nonzero norm. If the active edges
form a matching, two potentials can be counted by the
arbitrary-order one-grid estimate on an actual triangle.
The same tests apply to the individual least-repair
multipliers for \(q\mid M\), with a loss controlled by
\(1/(LM)\) independently of the number of edges. The
rank survives a common base change of sufficiently small
root discriminant, giving the full \(q^{3m^2}\) Kummer
extension degree. For any fixed \(k\) synchronized arrays
one can also keep their original edge observables and all
pairwise primitive secant odds in a single rank count.
With \(E_0=3m^2\), their joint Kummer degree is
\(q^{kE_0+\binom{k}{2}(E_0-1)}\) under the same field
condition. Taking odds cancels the normalization gcd exactly.
[Full power-class rank and least repairs (BX1–BX4)](ABC_SHARED_BOX_FULL_POWER_RANK.md)
· [Joint rank of synchronized arrays and primitive secants (SQ1–SQ3)](ABC_SYNCHRONIZED_PROJECTIVE_SECANT_RANK.md)

There is also an exact obstruction to universal polynomial
transfers. The only reduced three-term identities made of
edge-difference monomials on a simple graph come from a
triangle or the Plücker relation of a four-clique.
On the tripartite potential graph this leaves only the
original triangle relations. Nor can two distinct normalized
rational outputs both reuse only the original edge factors
if either output varies. On the actual boxes, a new bounded
degree irreducible factor generically contributes a
power-sized radical outside every original label. The
resultant bound is uniform over affordable catalogues of
polynomial-height coefficients. Several distinct new
factors contribute disjoint prime supports, with a
threshold \(t+s_i\) on a chosen active axis \(i\);
the original edge-private supports survive the same
forbidden-factor tests. This controls specialization as
well as universal identities.
[Polynomial transfers and fresh prime support (GP1–GP4)](ABC_SHARED_BOX_POLYNOMIAL_TRANSFERS.md)

## Private valuations, economical transfers, and real windows

For each fixed lower exponent on the omitted-vertex count,
every edge of a generic original box carries private full
prime powers with radical mass in a common compact interval
of positive scaled valuations. The interval and numerical
reserve are uniform in the packet and profile, and can be
chosen before the size of the box. On a fixed complete-letter
fibre, a prescribed multiplier \(Q\mid W\) also gives a
power saving for each interval of length \(QX^\theta\) in
the moving complete letter. These counts allow the interval
centres and multipliers to depend on the anchored values.
A forced covering by sufficiently few such intervals would
exclude that subfamily of boxes.
[Private valuation bands and cofactor dispersion (NB1–NB3)](ABC_SHARED_BOX_PRIVATE_VALUATION_BANDS.md)

For Laurent monomials in the original complete labels of
any fixed number of synchronized boxes, the private prime
powers bound the exponents after dividing by the actual
common gcd. For each residue list a coefficient dictionary
of size \(D\) then costs at most \(D/K\) of the spanning-tree
upper bound for nonidentity evaluations, up to a constant.
Thus an affordable dictionary cannot produce new primitive
additive outputs on a full-exponent cohort by those
monomials. The surviving identities give original
triangles or coefficient-only outputs. The bound on the total
exponent also allows slowly growing synchronized tuples and
dictionaries depending on inactive boxes. A \(K\)-sized
coefficient dictionary can make a nonidentity hold on
every box, so this threshold is sharp for the stated
primitive-height test.
[Laurent output count (LP1–LP4)](ABC_SHARED_BOX_LAURENT_OUTPUT_RIGIDITY.md)

The family derivative lattice has an exact spanning-tree
minor formula on private prime coordinates. Primitivity
makes the prime-depth components bipartite, so even very
large valuations contribute to the minor gcd only in
proportion to the cycle ranks of their depth graphs.
With the private numerical reserve above, the derivative
inequality has positive slack on sufficiently large fixed
boxes whenever the coordinates include a private prime
on each edge. This holds without a near-minimum packet
mass hypothesis and with arbitrary shared prime
coordinates added to the same matrix.
[Derivative minors and their height comparison (DM1–DM7)](ABC_SHARED_BOX_DERIVATIVE_MINOR_CYCLES.md)
· [Uniform application to substantial omitted-vertex exponents](ABC_SHARED_BOX_PRIVATE_VALUATION_BANDS.md#consequence-for-derivative-lattice-tests)

For the repaired target-support problem, a bounded-degree
polynomial carrier with small formal radical must have a
high-degree additive discrepancy. Keeping the contact
orders on the original formal support strengthens this
degree constraint. Slicing by integer lines turns the short
real window and exact local depth into a joint count,
uniformly in moving coefficients of polynomial height.
The reduction still needs control of actual target
multiplicities for thin families and for evaluations whose
radicals are much smaller than the formal degrees predict.
The Euclidean Fermat count for fixed signatures concerns
ordinary primitive rational images on fixed curves.
The uniform Fermat box estimates for three powerful
numbers start with primitive full letters. Neither
supplies the required least-repaired quotient
multiplicity estimate with varying supports.
[Real and local polynomial carrier counts (RW2–RW4)](ABC_POLYNOMIAL_REAL_WINDOW_OBSTRUCTION.md)
· [Achenjang–Arango-Piñeros–Vemulapalli, Theorem 1.1 and §5.1](https://arxiv.org/html/2611.31699v1)
· [Browning–Verzobio, Theorems 1.1 and 4.1](https://arxiv.org/html/2608.24512v1)

## Bounded occurrence tests and two-target reconstruction

An affine hash of \(a\) and \(c-b\) retains each true pair
with probability \(1/H\) and each false pair with probability
\(1/H^2\). Additional free buckets suppress false hits
without changing true retention. In a fixed complete fibre
the true pairs use distinct \(a\)'s, so pairwise independence
gives a constant-probability hit at multiplicity at least
\(H\). This replaces the multiplicity threshold by one
bounded occurrence test on the original eligible candidates.
With \(H=N^u\), additional bucket exponent \(\nu\), and
prime block exponent \(\beta\), a near-real candidate count
\(N^{1+o(1)}\) allows \(2u+\nu+\beta>1\) with a fixed margin.
The signed Jacobi version consequently uses a fixed ramp
threshold. A shorter numerical window permits a global
false-error bound by counting discrepancies only once.
[Hashed target occurrence (HT1–HT4)](ABC_HASHED_TARGET_SUPPORT_POPULARITY.md)

For a full equality grid, any test with \(d\) nonnegative
separated products and mean true acceptance \(p\) has
false rate at least \(p^2/d\) up to a finite-size correction.
Independent affine gates attain the scale. Requiring
all three original roles to separate gives the corresponding
\(p^{3/2}/\sqrt d\) scale on an addition table. On the
actual eligible incidence matrix, the largest \(d\)
multiplicities control the overlap with a rank-\(d\)
signed kernel. That energy bound ranges over the full
matrix on which the rank is given. A saving still needs a
bound restricted to actual eligible candidates.
[Separation rank and target energy (SR1–SR3)](ABC_RANDOMIZED_ADDITIVE_TEST_SEPARATION_RANK.md)

At fixed relative odds, common original divisors of the
two summands prescribe the quotient of two complete
targets. A mixed power defect and a numerical lock
of sufficient weight then recover the target roots by
projective spacing. This gives a strict count on diffuse
sectors where the target power defect is small, using
two target positions. For least repair, a residue of the
corrected multiplier ratio gives the analogous count at
any fixed positive numerical lock weight by reconstructing
the target radical ratio. This criterion places no
divisibility-grid restriction on the original valuations.
At minimum nonatomic mass, the corrected residue
dictionary must attain exponent at least \(1/6\) on a
saturating outside graph with locks of vanishing radical
price. The required uniform upper estimate on that
residue dictionary remains open.
[Two-target summand locks (TL1–TL5)](ABC_TWO_TARGET_SUMMAND_LOCKS.md)

A single original divisor \(D\) of positive numerical
weight also gives a direct reduction across the balanced
endpoint. Modulo \(D\), the corrected multiplier ratio for
the two other letters determines the \(M\)-th power of their
radical ratio. Since their radical heights are small, this
residue recovers the original row with subpower multiplicity.
The weighted number of residues **attained by actual rows**,
summed over a covering divisor family, has the same count
exponent as the primitive quotient family. Every balanced
endpoint has an allowable covering after choosing a positive
valuation band on a nonatomic letter. The new missing
estimate is a uniform power saving on those attained images.
[Direct primitive residue count (TL6)](ABC_TWO_TARGET_SUMMAND_LOCKS.md#a-direct-weighted-residue-count-for-primitive-quotients)

For targets with many pairs of representations, count each
normalized source pair once before restoring its actual
target completions. The resulting bound gives a strict
multiplicity-tail saving when the compressed shift price
and the source-pair codegree are both controlled. It keeps
the same target and all three original caps during
completion. At balanced individual cap \(2c_0\), with
\(c_0=1/6\), the sufficient inequality is
\(d+3v+3\rho<2c_0+3u\), where \(v\) bounds the codegree
exponent and each target provides at least \(n^{2u-\rho}\)
qualifying pairs of compressed price at most \(d\).
No uniformly cheap supply of such shifts is known.
[Compressed-shift target tails (SCM1–SCM5)](ABC_COMPRESSED_SHIFT_TARGET_TAILS.md)

## Shallow radical planes and deeper third-letter support

Three shallow corrected repair residues can force the
three residual radicals into a common rational plane
even when none of the individual complementary ratios can
be reconstructed. The congruence lattice controls the
full-rank occupancy. In a plane, the parts of the
three moduli absorbed by its primitive normal determine
the angular count and the cost of listing high-occupancy
normals. Actual primitive rows with squarefree residual
parts attain the absorption bound at ordinary radical
budgets. Thus primitivity and all three residues alone
allow the difficult normal content.
[Shallow records and the absorption example (SR1–SR4, SR6)](ABC_SHALLOW_REPAIR_RADICAL_PLANES.md)

On a radical plane, two equal exact residual powers
of order at least three give subpower completion even
with an arbitrary actual third valuation spectrum.
In the exceptional square case the third letter is
squarefree outside a fixed pool, and its positive
scaled-degree mass gives a strict saving on cheap
records. For unequal orders \(p>q\ge2\), the third
radical divides a binomial of degree \(p-q\) in the
second root. Fixing the exact quotient reduces the
original equation to one polynomial. Its factors
accounting for the third radical have bounded total
multiplicity, leaving a nonconstant factor supported
on a fixed prime pool.
[Equal-power completion (SR5)](ABC_SHALLOW_REPAIR_RADICAL_PLANES.md#two-equal-residual-powers-on-a-radical-plane)
· [Unequal powers (UP1–UP3)](ABC_UNEQUAL_POWER_RADICAL_PLANES.md)

For adjacent orders, high valuations on the original
third support give an exact quotient at any fixed
depth \(J\). With a suitable coefficient-height bound
\(n^\Theta\), second-root cap \(n^{r_y}\) and retained
third-radical floor \(n^{r_z}\), its normalized height
price is at most
\(\max(0,\Theta+r_y-Jr_z,\Theta-r_z)\).
Pay separately for low-valuation primes removed from
the third radical. The one-polynomial support
argument applies after this deeper quotient is fixed,
including for the small orders two and three.
The argument does not force either a radical plane
or economical records on general balanced profiles.
[Deeper third-letter quotient (UP4–UP5)](ABC_UNEQUAL_POWER_RADICAL_PLANES.md#a-deeper-quotient-forced-by-the-original-third-valuations)

For same-target pairs, the normalized source-pair
count now discounts the multiplier primes that remain
in the difference triple as well as those fixed by an
external pool. Once the two source values are known,
the shift exposes a prime pool for the actual complete
target and a second pool for the opposite summand.
Both lower the completion count through their
respective remaining radical caps. The combined
inequality gives a strict gain in parameter sectors
where the previously compared adaptive fixed-fibre
bounds reach their thresholds. It still requires
enough pairs with the indicated shift and overlap
prices.
[Retained support and three-letter codegrees (SCM8–SCM10)](ABC_COMPRESSED_SHIFT_TARGET_TAILS.md#endpoint-support-and-original-target-discounts)

## Divisible projected roots and simultaneous high-valuation tails

On a fixed projected radical plane, write the original
free letters as $a x^p u$ and $b y^p v$ with $u\mid x$
and $v\mid y$. For every valuation $e\ge p(p-1)$,
division of $e$ by $p$ gives this representation
without requiring a fixed valuation grid. The actual
third radical forces an integer quotient of the
corrections. A second integer separates the unmarked
root factors, and the original third valuations give
an exact quotient at any fixed depth. Once these
records are fixed, Mason–Stothers forces a nonconstant
polynomial factor supported on a fixed prime pool.
The same completion bound can charge one actual
correction in place of the separating integer.
[Divisible projected roots (DR1–DR5)](ABC_DIVISIBLE_POWER_ROOT_PLANES.md)

At unequal orders $p>q$, the smaller-order
correction divides the first quotient and can be
recovered at divisor cost. With affordable records
this yields a strict sector without a lower bound on
third-letter valuations. Multiple varying exact
layers supplied by a fixed polynomial carrier can
also be recovered from one quotient under the
given degree conditions. Both statements require
the indicated plane and record dictionaries.
[Unequal divisible roots (DU1–DU4)](ABC_UNEQUAL_DIVISIBLE_ROOT_PLANES.md)
· [Polynomial layers (LC1–LC3)](ABC_POLYNOMIAL_LAYER_RADICAL_PLANES.md)

In an actual completion star, a large least common
multiple of high-valuation divisors costs at most
$n^{2/R}$ to record, independently of the fixed
number of leaves. A degree bin with many neighbours
and large mean numerical tails consequently loses a
fixed power when the tails have small pair overlap.
Large overlap forces actual normalized shifts with
large common high-valuation divisors. The argument
retains the original target tail in the least common
multiple and bounds the admissible overlap through
the first two depth moments. A saving in the
remaining shift dictionary still requires additional
control of its radical.
[Simultaneous high-valuation tails (PL3–PL9)](ABC_POPULAR_COMPOSITE_PACKET_LCM.md)

## Affordable projected planes and differences after a common target lift

A strict saving follows on a family covered by sufficiently
few projected planes \(Ax+By+Cz=0\), where \(x,y\) are the
actual \(q\)-projected roots of two original letters and
\(z\) is the original third radical outside a recorded pool.
The saving permits \(q\lambda\to0\), every orientation and
arbitrary exact valuations. In the scale \(X^{1/q}\), the
third radical is a short multiplier. A cheap actual divisor
on one projected letter allows a positive-mass window on
the other to move, and original support recovery completes
the pair of letters. This removes the auxiliary quotient
prices from the plane route. A simultaneous residue branch
can supply the plane by a rank bound, with omitted-prime
multiplicity included when the branch is counted in original rows.
A single divisor also gives an unconditional dispersion test.
For the chosen window of original mass \(t\), a full-exponent
family needs at least \(n^{7t/8-o(1)}\) occupied mixed
divisor-residue records \((D,x/z\bmod D)\). Finding cheaper
dictionaries of that one corrected complementary residue
would suffice in place of extracting planes.
[Projected plane and mixed residue criteria (AP1–AP3)](ABC_AFFORDABLE_DIFFUSE_PROJECTED_PLANES.md)

For two primitive rows, lift to any common multiple of
their complete targets. The difference of the two lifted
opposite summands divided by their gcd is independent of
that common multiple. After cancelling the gcd of the
weighted cross products, the new difference is the reduced
shift times the original target gcd divided by the
original selected summand gcd. Thus fixed-degree norm
recovery charges a joint reduced-shift and power-class
record, plus the relevant original gcds. A common target
divisor or a short window in the relative odds can replace
the shift record by projective spacing. Both counts include
pairs with distinct least-repaired targets, and count the
primitive quotients without any common multiplier entropy.
Prime-index profiles show why limiting degree information
alone cannot produce a small joint defect.
[Primitive pair cancellation and counting (PD1–PD4)](ABC_PAIR_LIFT_PRIMITIVE_POWER_DIFFERENCES.md)
· [Fixed and moving target power-defect counts (NP1–NP10)](ABC_NORMALIZED_PAIR_POWER_DEFECTS.md)

Two independently moving targets can also be recovered
from a union of cubic defects without a prescribed exact
shift. A primitive diagonal cubic first recovers the
normalized opposite summands. A split cubic cover then
recovers their common rational scale. For three or more
rows with prescribed affine shape of the original
summand-to-target ratios, two split covers recover the
translation and scale from weighted joint power classes.
This replaces exact gap records by a shape record and
uses the actual valuations of all three letters.
[Affine shapes and joint power classes (AS1–AS6)](ABC_AFFINE_SHAPE_JOINT_POWER_CLASSES.md)

For actual completion stars, mix anchors from both summand
roles and charge their dictionary exponents separately.
The exact expected modulus is the weighted union of three
sources of prime powers, including the actual target tail.
A single anchor from the larger dictionary can supply weight
missing from many cheap anchors in the other role. The
cross moment also forces deep divisors on normalized cross
shifts. The remaining shift radical still needs an upper
bound. The two-moment threshold is sharp for abstract
occupancy data.
[Mixed-orientation completion stars (MSO1–MSO5)](ABC_MIXED_COMPLETION_STAR_OCCUPANCY.md)

## Small-modulus digit recovery and projective shape growth

Projected ratios at several orders give binary digit contrasts
of the original valuations. The contrasts are rational numbers
with signed squarefree support and survive a common unit twist
of the supplied residues. A logarithmic list recovers all but
vanishing radical mass at a modulus on the radical scale.
At any fixed smaller precision \(D>n^\delta\), a bounded
list of probes recovers nearly the entire nonatomic mass
on both free letters. The mixed Fourier estimate works with
atoms on either letter, since its resonant lines have zero
mass against the nonatomic carrier.

In a balanced box, one can choose an affordable full
divisor on a nonatomic opposite letter. From the recovered
diffuse supports, the remaining atom windows group into a
bounded number of roots and a short correction. The full
divisor then completes the original ratio at arbitrarily
small normalized cost. This supplies a bounded-length
small-modulus residue criterion for all balanced profiles.
The modulus and its full parent divisor can be recorded at
arbitrarily small fixed normalized cost. An argument forcing
a cheap *joint* dictionary for the resulting contrasts remains
necessary.
[Multiscale residues (MR1–MR4)](ABC_MULTI_SCALE_MIXED_RESIDUE_RECOVERY.md)
· [Binary support and atomic completion (MD1–MD3)](ABC_MODULAR_DIGIT_SUPPORT_AND_ATOMIC_COMPLETION.md)

For projective shape, a weighted product of actual secants
corrects the missing power class of the moving pole. The
corrected classes are invariant under separate common
multiplication of each row. Three rows require no shape record
in the cubic count. Small common original gcds reduce the
number of corrected classes needed to two. After fixing two
original anchors one can alternatively choose one uncorrected
class. Without such an anchor or correction, actual primitive
quadruples of fixed projective shape can share every even-weight
relative parity class in a power-sized family.
[Secant-corrected projective counts (PM1–PM8)](ABC_PROJECTIVE_SHAPE_SECANT_CLASS_COUNTS.md)
· [Two-class obstruction and uncorrected anchor count (PSO1–PSO2)](ABC_PROJECTIVE_SHAPE_TWO_CLASS_OBSTRUCTION.md)

Private secant primes on a fixed shape lie in the determinant
pool of its marks. After fixing one actual anchor, the private
packet on each of two secants puts the other endpoint into a
projective class within the available packing gap. Thus a fibre
of exponent \(r\) has at least \(n^{(s-1)r-o(1)}\) projective
shapes on its typical \(s\)-tuples, with the inner height
allowance chosen before the budget limit. This excludes the
cheap shape extraction on those typical affordable-fibre tuples,
even before charging the uncorrected power class. The estimate
does not extend to all specially selected tuple relations.

[Projective shape growth (PG1–PG3)](ABC_AFFORDABLE_FIBRE_PROJECTIVE_SHAPE_GROWTH.md)

## Positive residue products and moving powerful codegrees

A positive quotient of two projected powers captures the
same nonatomic carriers as the earlier parity contrasts.
The phase argument keeps a bounded list of masks below
any fixed radical-scale precision, even with atomic
partners. At each prime, two complementary ratios record
one of three exponent directions. In a two-row product
these directions cannot cancel. Exact reconstruction
therefore exposes the diffuse support of both rows, and
the affordable full divisors complete their atomic parts.
This gives essentially quadratic product image and
essentially diagonal multiplicative energy on balanced
charts. It also gives higher positive word bounds, with
subpower fibres at separated integer weights.
[Positive projection transcript products (PP1–PP4)](ABC_POSITIVE_PROJECTION_TRANSCRIPT_PRODUCTS.md)

The small-modulus congruences alone allow a much larger
contrast image than \(n^{1/3}\). A construction of actual
pairwise coprime integer rows satisfies the balanced
radical caps, simultaneous small additive congruences
on all three letters and a short real error bound.
Its contrast image can still exceed \(n^{1/3}\) by a
fixed power. The construction does not impose the exact
original sum. Thus obtaining enough multiplicative
collisions or an upper bound on the two-row product
image requires stronger additive information.
[Local additive transcript obstruction (SMO1–SMO3)](ABC_SMALL_ADDITIVE_MODULUS_TRANSCRIPT_OBSTRUCTION.md)

The omitted-coordinate graph of powerful sums allows moving
complete targets. A positive common target divisor recovers
a reduced complementary ratio, and a primitive diagonal
cubic then recovers the roots at subpower cost per
coefficient record. The target record charges the relative
cubic defect, with the extra shared support on a prescribed
integer. For the original row count, kernel multiplicity
enters only the linear term of the graph inequality.
In-column common target kernel factors give a strict
sector at the minimum dictionary sizes. They have no
global dictionary price.
[Moving-target powerful codegrees (MP1–MP6)](ABC_MOVING_TARGET_POWERFUL_CODEGREES.md)

On typical affordable-fibre triples the union of the
three corrected cubic records costs at least four times
the fibre exponent. Any two cost at least three times
that exponent, up to the fixed allowances. A separate
two-cap construction with diffuse exact-cube values has
unbounded normalized corrected cubic defects. Multiples
on the Fermat cubic also give power-many genuine primitive
tuples with both uncorrected classes trivial when the
original radical caps are absent. These bounds locate
which hypotheses matter in transferring the cubic counts.
[Corrected cubic class price and obstructions (CD1–CD4)](ABC_TWO_CAP_CUBIC_SECANT_DEFECT_OBSTRUCTION.md)

## All-row codegrees, aligned support, and exact products

The moving-target pair count extends to arbitrary moving
right values by dividing out their actual common gcd and
recording the outside relative power defect. The exact
determinant excludes repeated original rows and puts every
shared target prime in a fixed pool. Only the exclusive
target defect is paid. Consequently the all-row graph cost
can count pairs with equal left labels and drop the
cell-multiplicity term. The complete-label specialization
works at arbitrary original valuation profiles once its
exclusive target record is affordable.
[Moving-right and moving-target count (MT1–MT4)](ABC_MOVING_RIGHT_AND_TARGET_POWER_CLASS_CODEGREES.md)

For stars with several moving targets, the reduced root
ratios recover the original normalized coordinates. A single
weighted product class of the actual targets then gives a
split cyclic cover and determines the remaining common
scale. A divisor of the common omitted exponent suffices as
the class order. This counts all original rows with
repeated kernel labels allowed. The same three-cap
dictionary bound that controls anchored eligible pairs
also limits the gain from these stars. If the radical
masses at the required divisible valuations vanish, the
graph exponent reaches \(1/3\) even with favourable
spacing. Thus enlarging the stars alone does not recover a
fixed saving on arbitrary diffuse profiles.
[Joint target classes on actual stars (MS1–MS3)](ABC_MOVING_TARGET_POWER_CLASS_STARS.md)
· [Anchored extraction barrier (MK1–MK6)](ABC_MOVING_TARGET_KERNEL_EXTRACTION_BARRIER.md)

The exact products of the two complementary fractions
\(A/C\) and \(B/C\) over three rows have subpower fibres.
The nonsingular fibre is an elliptic curve with rational
\(3\)-torsion given by cycling the rows. At small external
modulus, positive quotient probes can separate scaled
valuation sums and avoid selected atomic windows. A prime
hidden from the three-row product then contributes only
a short correction to the exact products, or belongs to
bounded atom-membership groups. Two affordable actual full
divisors recover the two exact product ratios. This gives
essentially diagonal third multiplicative energy at
arbitrarily small radical-scale precision on the balanced
charts, at the cost of a logarithmically longer list of
projection orders.
[Exact complementary products (EC1)](ABC_EXACT_COMPLEMENTARY_PRODUCT_TRIPLE_COMPLETION.md)
· [Small-modulus and signed product completion (SC1–SC3)](ABC_SMALL_MODULUS_TRIPLE_PRODUCT_COMPLETION.md)

## Scalar product images and exact signed fibres

One scalar mask residue per probe, from the two free original
letters, gives arbitrarily small normalized completion cost
for positive words of every fixed length inside an affordable
full-divisor chart. A signed word admits the same completion
after fixing a set of rows that meets every zero-sum weight
subset. Alternatively its restriction to typical
small-overlap tuples has a small fibre without fixing an
entire row. Taking the positive word on tuples with a common
target expresses the weighted occurrence problem as a scalar
product-image count. For a target multiplicity threshold
\(n^u\) and a fixed \(k\)-fold positive product, a uniform
saving from \(n^{1/3+(k-1)u}\) for this image would suffice.
The ambient residue count leaves a power gap because the
direct lift pays at least the covered free-pair mass.
[Scalar products and target multiplicity (SR1–SR3)](ABC_SCALAR_MASK_TARGET_MOMENTS.md)

For the complementary projective record, finite interval
probes and Möbius inversion give exact signed word recovery
with the probe count bounded as \(L\to0\). With any three
distinct rows carrying two positive weights and one negative
weight, a prescribed exact word has fibre
\(O(n^{\rho/2+\eta}X^\varepsilon)\) under individual letter
radical caps \(n^\rho\). At most one role can share support
outside a fixed pool, and the two other caps give the saving.
This applies to primitive quotients in the original height
shell of a least-repaired count. The same signed completion
on arbitrary height-bounded rows has a rational-curve
counterexample.
[Finite probes and exact signed completion](ABC_FINITE_PROBE_AND_SIGNED_COMPLEMENTARY_COMPLETION.md)

Sequentially exposing stars in fixed right classes gives a
second constraint on cheap joint defects. Repeated support
or support already paid by a fixed group pool lowers the
remaining query caps. At vanishing aligned valuation masses,
this forces the existing star graph exponent strictly above
\(1/3\) for at least five rows if the free pool may include
left-label support and repeats, and for at least three rows
when its nonfixed primes must repeat on targets. The estimates
also apply in fixed complete-right groups. A pool assembled
from all pair determinants does not in general have the
required support property.
[Sequential support loss (MK7–MK8)](ABC_MOVING_TARGET_KERNEL_EXTRACTION_BARRIER.md#sequential-support-loss-in-larger-stars)

## Same-target word completion and cross-target residue counts

On a common complete target, the two exact complementary
products with mixed weights recover the target at divisor
cost. At zero total weight, equal products force a repeated
row when only one weight is negative. For \((1,1,-1)\)
the two denominators recover the target exactly and a split
quadratic recovers all three distinct rows at subpower cost.
For arbitrary fixed mixed three-slot weights, the free-role
support and the fixed-letter count give an exponent
\(\rho/4+\eta\) under individual letter caps \(n^\rho\).
The local positive complete-product image occupies a
subpower proportion of the allowed unit residue pairs
at all target prime-power depths. A count using only
those local congruences must allow this large image.
[Common-target complementary products (CT1–CT4)](ABC_SAME_TARGET_COMPLEMENTARY_PRODUCTS.md)

The subpower completion extends to all two-equal-positive
patterns \((h,h,-k)\) and to \((a,1,-1)\). The latter uses
the fixed input and output denominator in a homogeneous
resultant. For \((1,2,-3)\) with one complementary product
equal to one, an explicit elliptic curve with rational
three-torsion gives the same conclusion even for varying
targets. In the general balanced unequal pattern with this
exceptional output, the unresolved part of the binary
recovery is a common gcd on the opposite free role.
The exact word is counted at subpower cost once that gcd
is fixed, and its cost is at most its outside radical.
[Equal weights, rational maps and the exceptional gcd (CW1–CW4)](ABC_SAME_TARGET_EQUAL_WEIGHT_WORDS.md)

In a fixed target multiplicity bin \(n^{u+o(1)}\), one can
choose full-divisor lists shared by many complete targets.
The original witnesses force their residues into small
opposite-letter dictionaries. The summed exponent for
the number of complete targets is
\[
      E=\sum_i\bigl(\alpha_i-(u-t_i)_+
                       +\min(t_i,r_{W_i}-2u)\bigr)
\]
when the fixed list's least common multiple has sufficient
numerical weight. Here \(\alpha_i\) is the opposite complete
dictionary exponent, \(t_i\) is the packet's radical cost,
and \(r_{W_i}\) is the marked letter's raw radical cap.
Short real intervals reduce the necessary numerical weight.
A short-index progression estimate gives some further slack
when the step is large and has small radical. The actual
same-target scalar image has exponent at most \(E+ku+\eta\)
for a fixed \(k\)-fold word. This supplies a strict saving
in the bin when \(E<1/3-u\). The criterion includes intermediate
multiplicities and dictionaries larger than the minimum
balanced scale. Its compactness and numerical-weight
conditions still need to be forced where they are compatible.
In the tall balanced profiles, the same pricing would require
total packet radical mass below \(1/6\) for a strict exponent
saving. If the scaled orders available to all packets are
bounded by \(S<6\), their combined numerical weight is then
below \(S/6\), leaving a fixed gap from the strict CRT threshold.
[Full divisors and same-target image bounds (TC1–TC6)](ABC_TARGET_PACKET_CRT_MOMENTS.md)

For the scalar rational coordinate, exact \(r\)-th-power
targets alone allow \(X^{1/r-o(1)}\) distinct targets
with any fixed number of primitive square-summand
representations whose signed denominators retain radical
on the \(X^{1/r}\) scale. The Gaussian construction works
outside any prescribed polynomial-height prime pool and
requires no least-repair multiplier. It places no small
radical cap on the free letters.
[Exact-power target obstruction (GT1)](ABC_SAME_TARGET_GAUSSIAN_SIGNED_DENOMINATOR_OBSTRUCTION.md)

## Complete dictionary balance and joint residue cells

Fixing a cheap full divisor gives a square-root charge on
the number of permitted residue classes of an opposite
complete integer. At balanced one-third saturation this
forces the three complete-letter dictionaries to have
the same limiting exponent \(\alpha\in[1/6,1/3]\), and the
three complete degrees to have exponent \(u=1/3-\alpha\).
This holds across the balanced endpoint. On a cheap
fixed-divisor chart, typical residue cells have two
moving dictionaries of exponent \(u\) each and actual
row exponent \(2u\). A fixed dictionary imbalance gives
an exponent saving on charts with the required cheap
locks. A shared residue and height cell also reduces
the cost of repeated CRT witness dictionaries subject
to the packet and numerical height thresholds.
[Dictionary balance and residue cells (RB3–RB8)](ABC_BALANCED_DICTIONARY_RESIDUE_CELLS.md#equality-of-the-three-complete-dictionary-exponents)
· [Absolute residue charges (AR1–AR3)](ABC_ABSOLUTE_RESIDUE_TARGET_TAILS.md)
· [Shared witness counts (SW2–SW3)](ABC_SHARED_CRT_WITNESS_DICTIONARIES.md)

In same-target word arguments, the twofold coloured
collision count has a power saving over the abstract
product bound. The common-target threefold count also
saves on bins with multiplicity exponent \(u>5/36\), at
fixed distance from this threshold under the balanced
caps. At every depth large primes fill the natural joint
complete-product residue range when the free radical
caps are omitted. The missing input on regular residue
cells is a restriction on the two moving letters and
their actual signed sum or difference with all three
radical caps imposed.
[Coloured product energy](ABC_SAME_TARGET_COLOURED_PRODUCT_ENERGY.md)
· [Prime-power image](ABC_COMPLETE_WORD_PRIME_POWER_SURJECTIVITY.md)

## One-sided fixed-shift transfer into balanced cells

After fixing the actual joint deficient-factor record, a
single high valuation floor suffices for a strict
fixed-shift saving at its active height threshold. The
terminal cost is the joint dictionary of the remaining
low factors. It can be dropped if an independent full
divisor survives in the primitive shift, or if there is
a fixed relative real gap. In an ordinary residue cell
the marked divisor on both computed targets cancels
during normalization. Switching the common-row role and
the computed role retains the divisor, at a cost that
uses the actual marked-letter dictionary inside that cell.
Three simultaneous locks allow an extraction of positive
cells with all three complete dictionary exponents equal
on the local scale whenever the conditional row exponent
is positive. A local power saving below the derived
codegree threshold gives a strict row saving. All bounds
retain the original radical caps and repair eligibility.
[One-sided reverse transfer (RT1–RT6)](ABC_BALANCED_FIXED_SHIFT_REVERSE_TRANSFER.md)

Repeated copies of the same computed letter give
additional numerical tail constraints at every fixed
positive multiplicity. Outside the prime pool of the
anchor differences, high full factors can be packed
until their product exceeds the true query width.
These constraints and the optimized one-sided floor
test leave diffuse tall profiles. The compressed-shift
bounds require an effective price deficit on actual
same-target pairs. The residue identities do not
supply it, and generic reciprocal words on those
fibres retain large denominator radicals from their
private prime supports. Additive energy at the cubic
exponent still allows a subpower energy deficit, which
the existing qualitative exclusion does not rule out.
[Full-factor quantiles (CQ1–CQ4)](ABC_BALANCED_CELL_FULL_FACTOR_QUANTILES.md)
· [Additive, compressed-shift and reciprocal tests](ABC_BALANCED_FIXED_SHIFT_REVERSE_TESTS.md)

### Optimized chart prices and the surviving diffuse profile

A full-divisor chart that loses exponent \(f\) has a
codegree threshold \(u-2f\) on its full row bins.
Here \(u=1/3-\alpha\), and the anchors can be chosen
on the smaller free dictionary. A collision-weighted
version bounds the chart directly, with the actual
remaining radical on the anchor role and no cell size
mesh. The two-role record entropy charge also shows
that raw factor captures alone cannot win this codegree
comparison. Both estimates keep the actual chart loss
through the local count. The best comparison
still loses strength as \(u\) tends to zero.
For full-factor records of raw price, keeping the
original numerical height in the SHF floor test gives
a sharp variational obstruction. The wholly diffuse
tall \([2,4]\) profile has a strictly positive excess
above the local active height-floor threshold for
every \(0<t\le1/6\), even after optimizing free
stripping and paid excision. Raw half and nonatomic
tests alone optimize to \(\min(r/2,c)\) after
prescribed support stripping, with \(r,c\) the remaining
raw and nonatomic masses. An actual conditional
dictionary discount can improve those prices.

At shallow combined residue and real depth \(s\),
interpolation gives a full list-deficit charge on
fine residual layers when \(s\le c/2\), with \(c\)
now the remaining raw root cap. A joint record of
price \(b\) capturing raw mass \(f\) from one original
role gives a power saving if its conditional mixed
list exponent \(k\) satisfies \(s-k>b-f+u\) by a fixed
margin. This tests actual occupied labels even near
zero complete multiplicity. A joint interpolation
model attains the conditional information identities
with diffuse integer valuation dictionaries, without
imposing the additive and probe equations. The
additive-energy, compressed-shift and reciprocal
bounds leave open whether the simultaneous lifts
under all three caps can avoid a fixed list deficit.
[Exact chart price (RT9)](ABC_BALANCED_FIXED_SHIFT_REVERSE_TRANSFER.md#codegree-thresholds-at-the-actual-chart-price)
· [Collision-weighted codegrees and the raw-record barrier (RT10–RT11)](ABC_BALANCED_FIXED_SHIFT_REVERSE_TRANSFER.md#collision-weighted-codegrees-without-a-cell-size-mesh)
· [Sharp raw-floor envelope and joint probe obstruction](ABC_BALANCED_FIXED_SHIFT_REVERSE_TESTS.md#sharp-raw-active-floor-envelope)
· [Optimized occupied-list and collision bounds (JRM4)](ABC_JOINT_REAL_RESIDUE_MOMENTS.md#optimized-list-prices-after-full-factor-records)

## Induced matching and auxiliary tests at occurrence

At the complete occurrence endpoint, thinning a balanced
family gives an induced additive matching on the actual
complete dictionaries while retaining the exponent and all
three radical caps. A random prime product of normalized
size exponent \(2t\) preserves the exponent \(t\) of such
a matching. Recording a short real carry equation first
reduces the necessary modulus exponent to \(t\). The
remaining model uses a coarse integer coordinate as well
as the finite cyclic residue group. A uniform incidence
saving for this mixed relation would address the count
of occupied cells near zero multiplicity.
The real logarithmic moment also counts ordered pairs
of small-radical integers with no shared-support restriction
inside a prescribed relative strip. This separates the tall
labels and preserves the induced occurrence matching on
real bins alone. The spacing argument also works on the
balanced caps at an absolute height cutoff
\(X^{1-cL}\) for every fixed \(c>2/3\), and for \(c>1/3\)
when the tall dictionary used for spacing is wholly
nonatomic.
[Occurrence matchings and short real tests (IM1–IM3)](ABC_OCCURRENCE_INDUCED_ADDITIVE_MATCHINGS.md)
· [Real pair spacing (RSI1–RSI2a)](ABC_OCCURRENCE_REAL_PAIR_SPACING.md)

Two original radical caps and fixed full divisors on
the two capped letters allow full-exponent induced
arithmetic matchings, even with diffuse limiting free
profiles. One-letter polynomial evaluation dictionaries
can saturate their small-block root marginals and have
diagonal-scale product energy, while every affordable
full-divisor packet is determined by its complement.
A different false-candidate construction retains all
three radical caps, true least-repair heights and the
first two support-prime layers under additional
prescribed unit tests at arbitrary prime-power depth.
It has one dominant valuation degree and the polynomial
range for the auxiliary modulus has an explicit counting
cost. These constructions bound what the respective
combinatorial, projection and local tests can supply
by themselves.
[Two-cap induced matching (IM4)](ABC_OCCURRENCE_INDUCED_ADDITIVE_MATCHINGS.md#a-primitive-two-cap-induced-matching-with-diffuse-degrees)
· [Root-projection obstruction (RO1–RO3)](ABC_RESIDUE_OCCURRENCE_ROOT_PROJECTION_OBSTRUCTION.md)
· [Auxiliary unit tests (AU1)](ABC_OCCURRENCE_AUXILIARY_UNIT_TESTS.md)

## Arbitrary fixed coefficients and the reverse-transfer obstruction

A bounded collection of high roots can be released
with a separate integral cofactor on every letter.
A large complete value gives a coefficient large
enough for rational reconstruction of the two other
cofactors after the high roots are fixed. Combined
with radical-scale character moments and concentration,
this proves the two-cheapest bounded-group dictionary
bound with arbitrary prescribed coefficients. In
particular, a large primitive normalized codegree
quotient has the half-minimum bound without retaining
a unit shift modulus or paying for its low-degree
factors separately. Discrete compression to a bounded
collection of powered roots and an integral
\(n^{O(1)}\) cofactor gives subpower completion.
Paths with only small primitive quotients are
power-small at every fixed positive multiplicity.
Fix the first row and the gcd, and apply the
small-cofactor count directly to the second
original row with its common letter fixed. All original caps and
eligibility restrictions are retained.
[Boundary release and fixed shifts (MCB1–MCB3)](ABC_MULTIPLE_SHORT_COFACTOR_BOUNDARY.md)
· [Large and small primitive quotients (RT7–RT8)](ABC_BALANCED_FIXED_SHIFT_REVERSE_TRANSFER.md#arbitrary-weighted-shifts-at-radical-scale-height)

The full-range comparison uses the actual record
price \(b\), conditional free caps and numerical
heights. In a positive ordinary cell at fibre
exponent \(u=1/3-\alpha\), a conditional cap
reduction \(\Delta_i\) on one of the two free
roles saves if \(\Delta_i-2b>1/3-2u\), with
the normalized quotient height hypothesis.
Diffuse profiles on \([2,4]\) with saturated
small-block projection costs evade these tests at
every \(\alpha\). At zero multiplicity the missing
bound concerns occupied cells. A separate actual
record estimate applies at every multiplicity:
after fixing the three divisors it charges half the
sum of the two cheapest remaining quotient costs.
The same count bounds aggregate close pairs of
small-radical integers in all unit classes of an
arbitrary fixed step. In a single unit block of
fixed normalized index width the occupancy has
arbitrarily small normalized exponent, uniformly
in the step and the valuation profile.
[Full-range surviving profiles](ABC_BALANCED_FIXED_SHIFT_REVERSE_TESTS.md#surviving-conditional-profiles-across-the-full-range)
· [Progression pairs and fixed full divisors (PCP1–PCP4)](ABC_SHORT_PROGRESSION_PAIR_COUNTS.md)

A genuine common-neighbour path can be counted in either
original fixed-anchor fibre. The half-cap and nonatomic bounds
use the original heights, and the active-floor equality gives
compression without a low-degree terminal charge in the
derived shift. Fixed-rank sparsity makes outside quotients of
height at most \(X^\theta\) power-small for some fixed
\(\theta>0\) at each positive multiplicity threshold.
Conditioning the absolute residue bound on full-factor records
also forces their captured raw mass to be paid by information
about the original free value and its shared residue. This
condition still applies at zero multiplicity. Tall diffuse
profiles can pass the price tests even after arbitrary
determinant-support stripping at the allowed cost. The
additive energy, compressed-shift and reciprocal estimates
leave separate arithmetic inputs for those profiles.
[Original-fibre lift, conditional prices and surviving profiles](ABC_BALANCED_ORIGINAL_FIBRE_LIFT.md)

At the critical combined auxiliary test scale,
false induced matchings can have all three radical
caps and the exact target repair bound. Their
width is still larger than the fixed-power real
filter, and they fail the original exact addition.
A signed tree comparison yields polynomial local
mixed-character flatness for actual sparse residue
masks. A quantitative global lower comparison from
genuine triple histories is still required. In
bounded-degree polynomial models, an aggregate
depth criterion gives a power saving when the
resultant cost of the auxiliary list is controlled.
A growing list with uncontrolled resultant cost
can pass full-depth tests on false candidates.
[False mixed-bin matchings](ABC_OCCURRENCE_MIXED_BIN_FALSE_MATCHINGS.md)
· [Signed auxiliary comparisons](ABC_SIGNED_AUXILIARY_TREE_COMPARISON.md)
· [Aggregate support tests in polynomial models](ABC_POLYNOMIAL_AGGREGATE_SUPPORT_TESTS.md)

## Thin valuation windows and the occupied-cell problem

A gap among the finitely many retained orders extends
the arbitrary-shift boundary release to integral
cofactors of height \(X^{L^{1-\xi}}\) for a fixed
\(\xi>0\) chosen after the count allowance. At the
same scale one can retain groups with a growing list
of original exact orders in boundedly many thin
windows. Fixed-shift and arbitrary-step occupancy
then use these longer index ranges. The weighted
codegree argument at fixed positive multiplicity can
discard paths with primitive quotients in this
enlarged short range before applying the half bound
to the remaining window groups.
[Thin-window release (TW1–TW4)](ABC_THIN_WINDOW_ARBITRARY_SHIFT_RELEASE.md)

For a fixed record of three full divisors, slice the
actual complete-value lattice into parallel lines.
The largest admissible determinant deficit is on the
enlarged short index scale. Short-line occupancy counts
each line at arbitrarily small normalized cost,
including directions with uncontrolled radical.
For a scaled determinant \(\Delta\), the resulting
conditional bound is \(n^\eta(1+\Delta^{-1/2})\) up
to height-subpower losses. The scaled determinant uses
the two actual coordinate heights or the area of a
prescribed ratio window. This gives a saving without
a positive multiplicity assumption when the joint
record price and half its determinant deficit sum
to less than one third. A broad height-one diffuse
profile need not supply such a record. Using the lattice generated
by the actual rows also bounds their common projective-minor
gcd by the box area divided by the square of the row count, up
to count allowances. Empirical prime prices now charge only
half the numerical area deficit.
[Divisor-lattice lines and record tally (DL1–DL5)](ABC_DIVISOR_LATTICE_LINE_OCCUPANCY.md)

At occurrence, the complete dictionaries have a power
saving from the cubic additive energy exponent. A
fixed difference is split at a radical-scale primitive
quotient height. Small quotients place positive-scale
support inside their gcd. Large quotients admit the
arbitrary-coefficient half bound. Two of the three
marginals, and every mixed pair after role colouring,
have energy exponent at most \(5/6+\eta\). This gives
a quantitative \(L^1\) supply for centered auxiliary
masks. A finite-field induced matching can nevertheless
have ideal-scale marginal energies and character sums
and satisfy fixed power conditions, even at primes
above the quadratic occurrence scale. Nonnegative
separated tests at small complete multiplicity have
a degree-dependent false-mass lower bound on the full
product. These facts leave the signed joint or
arithmetic lifting problem.
[Complete occurrence energy and auxiliary masks (OE1–OE2)](ABC_OCCURRENCE_DICTIONARY_ADDITIVE_ENERGY.md)
· [Pseudorandom matching and separated tests (PI1–PI2)](ABC_OCCURRENCE_PSEUDORANDOM_INDUCED_MATCHINGS.md)

Projecting both auxiliary trees onto cycle actions
gives a product of the two cancellation factors.
For sparse masks and the stated overlap conditions,
this is a power below the unamplified cubic signal.
A direct construction that assigns one true row to
each split pins the child ratios to a small set and
is far from the required subgroup law. Thus a
remaining obligation is a joint arithmetic history
construction with widely distributed ratios and a
positive contribution from genuine rows.
[Two-sided cancellation and the direct split obstruction (CT1–CT3)](ABC_TWO_SIDED_AUXILIARY_TREE_CANCELLATION.md)

## Projective content, affine records, and auxiliary units

A small-radical primitive family has a quantitative
higher-order projective sieve bound. For a collection of
ordered distinct \(k\)-tuples with \(k\ge3\), let \(J^{k-2}\)
be the average number of extensions per occurring first
pair. The geometric mean of the common projective content
is at most the box area times \(J^{-2}\), up to the
height-first count allowances. Equivalently, the sum of
the weights of all prime-power projective coincidences,
averaged over the tuples, has a \(2\log J\) deduction
from the usual area bound. The same argument works in
angular windows anchored by the first pair. Empirical
frequencies of primes in actual projective minors price
divisors even on primes absent from the original letters.
At near saturation the price bound applies to the fresh
part on almost all ordered pairs.
[Projective content moments and the higher-order sieve (LC1–LC6)](ABC_PROJECTIVE_CONTENT_MOMENTS.md)
· [Empirical projective prices (PF1–PF6)](ABC_PROJECTIVE_FREQUENCY_DIVISOR_PRICING.md)

A prescribed affine line admits at most half the
least radical cost of its changing roles, with no
index-width restriction. Normalizing coordinate
gcds makes each two-role weighted difference primitive.
This yields a second successive-minimum test and an
inverse-square bound for the triangle-minor gcd in
translated boxes above the affine-line floor. The
divisor-lattice dual has explicit frequencies in
terms of the three divisors. Raw prices in the
height-one \([2,4]\) diffuse example are still too
large to afford a dual frequency of small range.
Separately, two compressed complete letters can
now be released together without fixing the third
one. In particular two divisors leaving short
quotients give a subpower count per joint record at
every complete multiplicity.
[Affine lines and translated content (LA1–LA4)](ABC_LONG_AFFINE_LINE_RADICAL_COUNTS.md)
· [Two compressed letters (DTC1–DTC2)](ABC_TWO_COMPRESSED_LETTER_RELEASE.md)

A binary-form projective pullback bound shows that
the false Pythagorean carriers can pass the first two
support layers and still avoid common projective
clustering at every fixed positive numerical scale.
For a different route, the character large sieve
supplies nearly uniform independent spectator units
from the actual one-role radical dictionaries.
A bilinear transformation inserts such auxiliary
parameters at a single true-tree split, while
multiplying the exact additive defect by just one
parameter. Its output marginal depends only on the
normalized defect. The signed argument still needs
a way to retain the original labels at low cost
and to coordinate this insertion in two histories.
[Projective pullbacks and false carriers (PP1–PP4)](ABC_PROJECTIVE_PULLBACK_CLUSTERING.md)
· [Multiplicative unit supply (MSU1–MSU2)](ABC_MULTIPLICATIVE_SPECTATOR_UNIT_SUPPLY.md)
· [Defect-preserving auxiliary splits (BS1–BS2)](ABC_AUXILIARY_BILINEAR_SPLIT.md)

## Shallow occupied registers and short real auxiliary scales

A residue-list moment can include collapsed root products
larger than the modulus, with their actual multiplicities
in a congruence class. Interpolating the moment orders shows
that, at zero complete multiplicity, any joint unit probe of
depth at most \(1/6\) on one free letter has full collision
entropy, even with macroscopic exact valuation layers.
When the individual layer caps tend to zero it reaches
the full complete exponent \(1/3\) at depth \(1/3\).
Truncate three genuine locks to depth \(1/6\) each. In
that fine-layer regime, three shallow residue alphabets
then carry an actual partial Latin square with exponent
\(1/3\) even at zero complete degree. Prime-power probes
also give multiscale nonconcentration on a retained
marginal. The lifting problem now asks which dense
shallow arrays can have representatives satisfying
the three original height and radical constraints.
A CRT example lifts a complete Latin square to
primitive tall rows with the prescribed locks when
the complete radical caps are removed.
[Interpolated shallow residue moments](ABC_SHALLOW_MOMENT_INTERPOLATION.md)
· [Shallow registers and actual rows](ABC_SHALLOW_RESIDUE_REGISTER_SATURATION.md)

For the signed auxiliary comparison, the weighted
character frame applies to a fixed product of auxiliary
primes. A Pell ratio close to one changes its selected
quadratic sign and preserves the real-window bound to
power accuracy. A ratio of independent products supplies
central-window mass by pairing products in the same
short log interval. On paired blocks of signed histories
this provides actual short real leaf scales and retains
the exact paired products, with original-label recovery
from the normalized split data. The signed upper bound
now holds with those real windows. A positive signed
lower bound from genuine three-letter incidence is still
required, and independent spectator replication loses
too much of that signal.
[Simultaneous CRT frames and Pell correction](ABC_SIMULTANEOUS_CRT_CHARACTER_FRAMES.md)
· [Original-label recovery in the signed comparison](ABC_SHARED_SCALE_BILINEAR_TREE_HISTORIES.md)

## Joint real bins, numerical splices, and prime interpolation

A coupled list of unit residues and real logarithmic bins
has the same moment gain as one probe with their depths
added. The proof keeps the spacing of collapsed integer
root products in each torsion residue and uses a finite
Fejér majorant for the logarithmic bins. On a fine-layer
zero-multiplicity chart the mixed collision entropy
saturates up to depth \(1/3\). After recording a full
factor it saturates up to the remaining raw mass. This
rules out, by a power saving, actual charts on which a
residue label and the fixed records confine enough of a
complete role to too few numerical windows.
[Joint real and congruence moments (JRM1–JRM3)](ABC_JOINT_REAL_RESIDUE_MOMENTS.md)

In a shallow register of depth \(\sigma<d\), a matching
of exponent \(d=1/3\) produces \(n^{2d-\sigma-o(1)}\)
compatible cross-row sums or signed differences. Their
outside output radical generically exceeds the original
target cap. This remains true for numerical near-splices
at the mixed precision, with the precise decrease in the
radical threshold given by the real depth. The complementary
CRT examples illustrate two limits: a near-full shallow graph
can have no spurious triangles, and a full shallow Latin
square can satisfy the mixed entropy conditions while its
spurious carries miss exact addition by the bin spacing.
These lifts retain the full lock and primitive conditions
and omit the complete radical caps.
[Splice radical threshold and carry examples (SRS1–SRS2)](ABC_SHALLOW_REGISTER_SPLICES_AND_CARRIES.md)

Random interpolation through actual triples, combined with
the polylogarithmic-degree curve saving, rules out heavy
joint residue compression over a large logarithmic weight
of prime tests. The subset of rows exhibiting compression
can vary at every prime. The same argument bounds a proper
joint collision statistic over prime powers. The marked
lock primes alone supply too little weight to force the
obstruction. In the \(p\)-adic projection direction, the
original complete-ordinate graph at sparse terminal depth
is below the density threshold even under optimal marginal
power estimates. A transfer would require new variables
or a dense restricted graph with a controlled image.
[Joint interpolation and projection limits (JSI1–JSI4)](ABC_SHALLOW_JOINT_INTERPOLATION.md)

## Common secant gcds and private supports in hybrid cells

The joint real–residue moment gives a one-class saving
equal to the combined depth at low radical-scale depths.
In a regular near-occurrence chart, this bounds the full
part of a common secant gcd on a prescribed prime pool.
Adaptive pools from earlier anchors can be combined
with the interpolation hit bound. At zero complete
multiplicity the resulting mean normalized log gcd
tends to zero, and in sufficiently large shallow cells
the imposed modulus accounts for the mean gcd up to
vanishing extra depth. An ordered-anchor argument also
gives a fixed-budget power tail above each prescribed
positive normalized depth. Its power loss may shrink
with the budget. Separately, the part on prescribed
pools and the excess powers at small outside primes
admit a uniform power tail at fixed positive excess.
[Mean secant gcd (SGM1–SGM4)](ABC_SECANT_GCD_MEAN_OVERLAPS.md)
· [Normalized tails and primitive directions (NGT1–NGT3)](ABC_SECANT_NORMALIZED_GCD_TAIL.md)
· [Uniform tail for powerful excess (SGP1)](ABC_SECANT_GCD_POWERFUL_EXCESS.md)

With small scaling gcd, a fixed primitive direction
fixes an original complete-value shift up to a small
scalar catalogue. The fixed-target half bound controls
its fibre. This prices directions below an explicit
total radical threshold of \(3/2\) near occurrence,
apart from the fixed-budget exceptional set. Stronger
qualitative one-coordinate exclusions are available
from the unrecorded-cofactor method. The secant gcd
itself therefore need not account for the large radical
of a typical primitive direction.

On an actual half-depth hybrid cell a prescribed pool
containing an appreciable part of a free letter's radical
also lowers its per-anchor count. This makes the two
opposite supports almost private within genuine large
shared-residue cells in the stated fine-layer range and
permits recovery from reciprocal denominators. For a
numerically near splice, its gcd with the compared true
target has only a small excess over the shallow
modulus below half depth. Up to full depth, the averaged
factor test still excludes large supported or powerful
excess. The squarefree outside part carries the full
record price. The resulting direction and splice
estimates describe the obstruction to transferring the
original family into a smaller count of powerful
primitive quotients.
[Private support and splice gcds (HPS1–HPS3)](ABC_HYBRID_CELL_PRIVATE_SUPPORT_AND_SPLICE_GCD.md)

## Exact common factors and typical physical shifts

In a nonprimitive derived relation, primes of the common
gcd that survive in the primitive coordinates are already
visible there. The missing primes require a separate
record. Their exact radical credit pays for this record
at the one-third coefficient. The same compensation works
for two charged roles. For several additive branches
sharing one original value the union of primitive supports
restores the common value only once, leaving a smaller
missing-support charge.
[Exhausted-support restoration (ESL1–ESL5)](ABC_SPLICE_EXHAUSTED_SUPPORT_LIFT.md)
· [Common-anchor branch counts (SAB1–SAB2)](ABC_SHARED_ANCHOR_PRIMITIVE_BRANCH_COUNTS.md)

Near zero complete multiplicity, the mixed moment bound
controls the **pair mass** of oversized cells at every
fixed subterminal depth. Large exact valuation layers
can be split at their original exponents when individual
primes are small. Recording the large individual primes
preserves the nonatomic valuation mass in a full residual
chart. On the resulting charts same-role overlaps of
fixed positive mass are power-rare. The physical shift
and the actual positive splice then each escape every fixed
outside-radical budget in the stated pair populations.
For the normalized comparison with the marked target,
a positive average extra gcd is also power-rare across
full cells. These statements are compatible with large
actual fresh prime supports, which a closing transfer
would have to control by a different arithmetic input.
[Fine subroots and pair pruning (JRM5–JRM6)](ABC_JOINT_REAL_RESIDUE_MOMENTS.md#fine-residual-caps-from-small-individual-primes)
· [Generic large shifts and splices (HPO1–HPO4)](ABC_HYBRID_SAME_ROLE_OVERLAP_AND_FRESH_RADICALS.md)
· [Recurring splice divisors (RSC1–RSC2)](ABC_RECURRING_SPLICE_GCD_CELL_COST.md)

The optimized mandatory repair for the formal broad
\([2,4]\) profile leaves a quantitative scalar gap to
bridge, even allowing unequal fullness orders. The
adaptive valuation basis has an independent tradeoff
between the number and sizes of moving root coordinates
and the numerical width of a correction.
[Repair and valuation basis tests (DB1–DB3)](ABC_DIFFUSE_REPAIR_AND_BASIS_OBSTRUCTIONS.md)

The short-power-window construction gives primitive
two-cap triples at the limiting \(1/3\) exponent,
even with several pairwise coprime near-power occupants
per centre. It also permits a growing full prime-power
divisor of the displacement at fixed fullness order.
Its numerical lock weight vanishes when the budget
tends to zero, and it leaves the whole displacement
radical unrestricted.
[Coprime short-power windows (CPW1–CPW2)](ABC_COPRIME_SHORT_POWER_WINDOWS.md)

A quantitative Subspace Theorem application gives an
unbounded-height cardinality estimate for weighted powers
with a fixed finite set of places. The finite-height shared-order sector uses a separate
elementary count. After fixing the coefficient and third-divisor
records the addition fixes a unit power residue of the primitive
slope. Positive cross products in one class coincide as soon as
\(D\ge X^{2/k}\). This gives height-subpower completion uniformly
over all adaptive orders and primes dividing them. In particular a
divisor of fixed height \(D\ge X^w\) allows every threshold
\(K(L)\to\infty\) at late fixed budgets, and fixed orders with
\(kw\ge2\) under exact height bounds. Large orders are counted by
enumerating the small roots. The actual joint record price controls
the sector saving.
[Shared-order sector (CPW5)](ABC_COPRIME_SHORT_POWER_WINDOWS.md#an-adaptive-common-power-sector-with-a-third-letter-lock)
· [Unbounded-height statement (CPW3)](ABC_COPRIME_SHORT_POWER_WINDOWS.md#uniform-upper-bounds-once-the-shift-part-is-larger)

## Primitive slopes and radical prices of power differences

There is a finite-height improvement over the common-order
application above. The actual two powered letters may have
different literal orders. If their orders have a common divisor
\(g\), count the resulting primitive slope in a unit class
of a fixed full third divisor. In a true difference window
the count charges the actual real width and the product
of the two remaining root heights after stripping known
factors. A fixed \(g\) with \(gw>2\) suffices for subpower
completion at third-divisor height \(X^w\). This removes
the growing-order requirement for the finite-height application
with the same strict joint record price.
[Primitive slope count (CPW7)](ABC_COPRIME_SHORT_POWER_WINDOWS.md#finite-height-slope-counts-for-unequal-orders)

APC2 balances the original three-letter prime-budget count
against reconstruction from a light projected coefficient record
and the actual third divisor. At divisor price \(\mu<1/3\) its
optimal split is \(t=(1-3\mu)/(1+2a)\), giving exponent
\(b_0+[a+(1-a)\mu]/(1+2a)\). In particular \(a=1/2\) gives
\(1/4+b_0+\mu/4\). APC7 lets a prescribed sparse prime pool
also supply part of the third divisor. The same formula uses
\(u=\mu_0+\beta\), where \(\mu_0\) is the outside divisor
radical price and \(\beta\) is the actual conditional ratio
dictionary price after fixing the divisor and projection.
In centered fine grids the coefficient is recovered from the
actual root ratio modulo the divisor when \(D\ge X^{kL}\).
Thus the conditional record can instead describe the root ratio,
including a sufficiently small relative-height cover.
[Optimized coefficient pools (APC2)](ABC_ADAPTIVE_POWER_COEFFICIENT_PRICES.md#a-separating-pool-for-the-joint-coefficient)
· [Conditional joint pools (APC7)](ABC_ADAPTIVE_POWER_COEFFICIENT_PRICES.md#conditional-separation-with-a-joint-divisor-pool)
· [Root-ratio recovery (APC8)](ABC_ADAPTIVE_POWER_COEFFICIENT_PRICES.md#conditional-recovery-from-the-root-ratio)

A different record keeps only the mixed fourth-power defect of
\(UVW^2\). Its support detects residues modulo four on the
two operands and odd valuations on the third letter. Fixing the
defect gives a uniform bounded-height count on an elliptic curve
with rational two-torsion. A sparse pool can price this joint
defect at the same optimized exponent with its outside defect
price in place of the divisor price. No numerical lock is needed.
A saturating family must have mixed defect price at least
\(1/3-o(1)\) in every orientation, but no cheap conditional
defect record is forced on enough rows.
[Mixed fourth-power defect records (APC9)](ABC_ADAPTIVE_POWER_COEFFICIENT_PRICES.md#a-mixed-fourth-power-defect-replaces-exact-power-records)

An independent curve argument allows the centered signed
coefficient ratio to move without recording its dictionary.
It needs small numerical coefficient height and a third divisor
with affordable support and controlled scaled valuations.
The full powers give enough congruence depth to force a
plane equation. High-degree components use a uniform
rational-point bound. Low-degree components have few possible
roots modulo a sufficiently large part of the divisor. A
fine-grid averaging calculation quantifies the numerical
coefficient price on a diffuse spectrum. It does not force
the coefficient width needed at the same order and modulus.
In the separating-pool count, the pool need only lead
to few small relative-height classes after the actual
third divisor and projected valuation record are fixed.
A weighted monomial count then handles the coefficients
inside each class. If \(\mu\) is the divisor price and
\(b_1\) the conditional cover price, put \(u=\mu+b_1\).
The optimized exponent is \(b_0+[a+(1-a)u]/(1+2a)\).
For high-depth divisors this preserves the separating
exponent when \(b_1=0\), and also works with comparable
unequal orders by regrouping their exponents.
[Conditional relative-height covers (APC4)](ABC_ADAPTIVE_POWER_COEFFICIENT_PRICES.md#price-of-a-cover-by-small-relative-coefficient-heights)
· [Moving coefficient count and the centered benchmark (MPC1–MPC2)](ABC_MOVING_SIGNED_POWER_COEFFICIENTS.md)

A relative-height cover in a true near-saturating
family must respect arithmetic support overlap.
Remove atomic slots from two stabilized roles and
let \(M\) be their combined nonatomic mass. At
typical fine-grid orders the retained centered
coefficients have generic relative height exponent
\(M/4\) in units of \(kL\). For a fixed reference,
neighborhoods of exponent \(c<M/4\) cost at most
\((1-M+2\sqrt{cM})/3\) in base \(n\), up to a count
allowance. The two most diffuse roles have \(M\ge1/3\)
at every balanced limiting profile. The arithmetic
matching construction with two capped operands
still needs the full coefficient cover price at
small radii. A proof of the unrestricted saving
requires a conditional extraction that beats these
constraints and loses less in retained rows than it
saves in record price.
[Fine-grid relative distances (APC5–APC6)](ABC_ADAPTIVE_POWER_COEFFICIENT_PRICES.md#fine-grid-distances-at-balanced-saturation)
· [Numerical cover obstruction](ABC_ADAPTIVE_POWER_COEFFICIENT_PRICES.md#numerical-relative-height-covers-on-the-same-matching)
· [Quantified three-cap extraction obligation](ABC_ADAPTIVE_POWER_COEFFICIENT_PRICES.md#the-extraction-obligation-with-all-three-original-caps)

Each nonlinear rational factor of a fixed-order weighted
binomial pays for the original primitive roots at its own
outside radical price. Linear factors give genuine additional
additive rows with controlled radical caps. Averaging these
prices saves when the binomial has sufficiently many rational
factors and the coefficient record is affordable. For differences
of exact powers of order \(M=2^m3^b\), the global sector
coefficient is at most \(1/(\tau(M)+2)\). A separate bound uses
only the original shift cap in the numerator and \(\tau(M)\)
in the denominator. Full divisors of the shift are optional.
[Nonlinear factor and weighted-binomial counts (BFR1–BFR2)](ABC_BINOMIAL_FACTOR_RADICAL_PRICES.md#individual-nonlinear-factor-fibres)
· [Exact power-difference compression (RC1–RC5)](ABC_EXACT_POWER_DIFFERENCE_RADICAL_COMPRESSION.md#statement-and-counting-input)

Fresh cyclotomic factors can also be evaluated on any original
primitive row. In one orientation a fixed collection of \(J\)
such outputs with union radical at most \(n^t\) outside the
original triple and a prescribed pool gives a count coefficient
at most \((1+t)/(3+J)\), up to an arbitrarily small allowance
at small fixed budgets. The simultaneous version uses outputs
from all three orientations and counts their common order-six
factor once. These inequalities force new radical growth in
any near-saturating family across the entire dictionary range.
An upper bound for that fresh radical on enough genuine rows
would give a saving without a valuation-grid hypothesis.
[Fresh cyclotomic radical counts (BFR3–BFR4)](ABC_BINOMIAL_FACTOR_RADICAL_PRICES.md#fresh-cyclotomic-radicals-on-original-rows)

## Fixed-letter pairs and independently supported anchor stars

When two rows share a fixed complete letter, the
relative quotient of their complete free odds is
either in a deep unit congruence or in a small real
window. A two-base count bounds the exact quotient
when the full power coefficient ratios are numerically
close. The split cross-ratio identity then recovers
both rows at subpower cost. This gives a subpower
bound on all close pairs per fixed letter, including
short fixed summands. The reduced quotient also recovers
a moving common target at divisor cost. A prescribed deep
divisor shared across targets therefore replaces the whole
complete-target dictionary in the count. In an actual switched residue
cell the same recovery applies to pairs of queries
after fixing two anchors. The chart divisor and the
anchor difference together give a modulus containing
the square of the original divisor. The normalized
operand gcds can supply additional usable powers.
[Fixed-letter and switched-cell pair recovery (FRP1–FRP5)](ABC_FIXED_FIBRE_RELATIVE_ODDS_PAIRS.md)

A complementary bound uses centered coefficients on
selected diffuse slots on any fixed list of complete
anchors. After stripping the primes of anchor
differences, the computed values have disjoint
supports at each query. Short relative-height pairs
then constrain the residual radical mass in *both*
directions on each original fibre. With \(h\)
selected roles of residual raw cap at most \(R\),
the power in the per-query pair count is at most
\(2\sqrt{Rc/h}\) for diameter exponent \(c\) in
units of \(kL\). The separate raw and nonatomic
caps give the sharper profile-dependent formula.
This establishes a codegree saving criterion at
every fixed positive multiplicity when the actual
conditional cover price is sufficiently small.
The fine-grid observable can omit atomic slots in
this count. Exact two-base recovery instead needs
the entire quotient or a separate bound absorbing
the omitted factors.
[Common-query star quantiles (ACQ1–ACQ4)](ABC_ANCHOR_STAR_COEFFICIENT_QUANTILES.md)

Exact common powers have a uniform strict saving
for every adaptive common order, on either side
of the equation. The same-side quadratic estimate
covers orders that are powers of two. Capelli's criterion reduces any rationally reducible
weighted pure binomial on two primitive letters to
those sectors or the fixed Sophie Germain form
\(x^4+4y^4\). This gives a strict saving with arbitrary
coefficient dictionaries on the entire reducible
sector. A separate two-cap construction shows that
every fixed collection of fresh cyclotomic outputs
can have large radical even on an actual diffuse
complete-letter matching. Its third original
radical violates the small-budget cap by a growing
margin.
[Adaptive and reducible binomial sectors](ABC_BINOMIAL_FACTOR_RADICAL_PRICES.md#exact-power-sectors-with-arbitrary-common-order)
· [Fresh-output obstruction](ABC_BINOMIAL_FACTOR_RADICAL_PRICES.md#a-two-cap-obstruction-to-any-fixed-fresh-output-catalogue)

## Complementary words and reduced least-repaired target ratios

In the occurrence sector, two fixed positive-depth divisors on
different original roles give two congruences for the complementary
quotients of any pair. If both quotients have full narrow-root and
short relative coefficient representations, the exact outputs
recover the two rows even when their complete targets differ.
The count on that chart is \(O(n^\eta X^\varepsilon)\) for any
fixed \(\eta>0\) under the stated width conditions. One divisor
also gives a saving for close pairs whose free scaled valuations
are bounded, because the nontrivial reduced odds quotient must
carry a fixed amount of radical. Both arguments impose all three
original caps. The missing extraction is a power-large supply of
jointly close pairs on affordable original two-divisor charts.
[Cross-target complementary quotients (CTO1–CTO2)](ABC_CROSS_TARGET_COMPLEMENTARY_WORDS.md#two-complementary-locks-recover-all-eligible-pairs)

At every fixed order, exact two-coordinate signed product fibres
have a strict saving over their half-length radical price.
Outside the prescribed target support, primes occurring on just
two rows must occur on the same role with opposite signs.
Reserving the largest packet for exact completion gives the
saving. For equal positive weights the first possible packet
uses three rows and the three-row elliptic count strengthens
the bound.
[Signed complementary word counts (SWP1–SWP2)](ABC_CROSS_TARGET_COMPLEMENTARY_WORDS.md#exact-signed-complementary-words-with-all-three-caps)

One exact primitive odds quotient and the reduced ratio of the
two target radicals recover the pair at height-subpower cost.
A norm record on the target radicals or a short numerical
window on the actual least-repaired ratio can give the second
record. There is an obstruction to using cheap whole-difference
records on typical pairs. If a family of exponent \(d\) has
per-target-radical multiplicity exponent at most \(u<d\), a
repair order divisible by a sufficiently composite fixed integer
forces arbitrarily large outside radical on the normalized
repair differences of almost all neighbors of each reference
row. The argument is uniform over prescribed polynomial-height
reference pools. For actual least repairs it suffices to fix a
common wild tail, with a cost that can be made arbitrarily small.
The obstruction holds throughout the balanced multiplicity range.
[Repair ratio recovery and fresh radical (RRP1–RRP4)](ABC_REDUCED_TARGET_RADICAL_ODDS.md)

## Expansion of occupied residues and joint ray counts

For a deep absolute modulus, an occupied residue set \(S\)
can be charged by the square root of \(|S\Theta|/|\Theta|\)
for any prescribed auxiliary set with product of size at most
a fixed power of \(n\). The relative expansion majorant fits
the same amplified character moments as the absolute residue
count. On a balanced chart of dictionary exponent \(\alpha\)
this forces near-maximal expansion against every such auxiliary
set for each full-size occupied list. All fixed-length signed
modular words of those lists therefore have full image exponent
\(m(2\alpha-1/3)\). Joint valuation slots also bound every fixed
product moment at the diagonal exponent and bound the image of
arbitrary selected tuples. This obstructs residue compression
throughout the positive-surplus range,
including occurrence. The minimum dictionary has surplus zero.
[Relative expansion (RX1–RX2)](ABC_RESIDUE_EXPANSION_AND_JOINT_RAY_COUNTS.md#every-large-balanced-residue-list-has-maximal-finite-order-growth)
· [Diagonal modular energy (RX7)](ABC_RESIDUE_EXPANSION_AND_JOINT_RAY_COUNTS.md#signed-residue-words-and-diagonal-energy-on-the-cells)

For prescribed quadratic parity fields, a common target lock
also couples occupied residues of several actual rows in a
joint ray quotient. Each target prime of equal positive
order is priced once across the rows, and its radical price
is halved. The comparison of short ideal generators uses
the two different quadratic fields in their compositum.
Both the normalized ratio and the actual absolute residue
versions retain all original caps. In a composite packet
block the resulting gain can come from a smaller joint
field catalogue, from equal-order target overlap, or from
target primes in another row's field kernel. The bound
charges the actual price of finding those pairs, and
at the independent field count and saturated raw
caps, the cost of extracting any single-direction
overlap exhausts that ray saving. The estimates
then require overlaps in both directions or an
improvement in the actual joint field price.
[Joint target ray count and block criterion (RX3–RX5)](ABC_RESIDUE_EXPANSION_AND_JOINT_RAY_COUNTS.md#joint-target-ray-groups-and-shared-radical-prices)

## Joint normalized modular records and cross phases

Two products of original letters give the monic quadratic of
the normalized residues on each full-divisor lock. The resultant
of two matching such polynomials amplifies the modulus dividing the
cross-determinant product. At precision two, total numerical lock
weight strictly above one recovers a complete unordered pair.
At any fixed positive lock weight, sufficiently fine precision
also recovers exact complementary products. The normalized word
sieve can rotate through the three cap orientations and charge
prime overlap across different original roles. A tight product
of full locks also limits the actual integer shift dictionary,
giving a strict sector at positive multiplicity when the third
residue register has sufficiently small additive difference set.
[Normalized records and locked shifts (NJ1–NJ5)](ABC_NORMALIZED_JOINT_MODULAR_MOMENTS.md)

For a pair in distinct parity fields, an integral cross phase
uses only their symmetric-difference squareclass. A real quadratic
ray count charges the common target radical and one reduced
shift at half price, with an explicit extra charge on repeated
operand primes. A quotient of the two cross phases cancels the
common target factor and the unknown absolute residue, leaving
a half charge on the two reduced shift radicals and their
exceptional overlaps. At a common parity field and an equal
odd target, the two cross phases instead give two genuine
primitive abc rows. Their square targets partition the original
target prime powers. Their operand powers come from the two
actual reduced shifts. This transfers both exceptional-set
and high-fullness estimates to the pair when those shift
records are affordable in the original family.
[Cross-phase pair counts and auxiliary triples (CP1–CP5)](ABC_CROSS_PHASE_AND_AUXILIARY_TRIPLES.md)

## Joint displacements and one-auxiliary recovery

For a fixed vector of exact changes to the three complete
letters, each role gives its own fixed-shift dictionary.
Pruning the translated chart and combining all three residue
supports bounds the joint fibre below the two-letter count
in a range of positive multiplicities. On full-scale genuine
pairs, small primitive quotients can be discarded by the
small-feature bound, so the remaining shifted projection
prices are half their individual outside caps without a
low-valuation restriction. The count prices the actual
joint signed residue differences or an independent exact
shift catalogue. At a deep three-lock chart the pairs
with all three registers equal alone exclude the
maximal-depth boundary throughout \(1/6<\alpha<2/9\),
with no additive-compression assumption.
[Displacement fibres and the deep same-register bound (JD1–JD4)](ABC_JOINT_DISPLACEMENT_CATALOGUES.md)

The determinant of the products of the two scaled operand
pairs gives one ordinary primitive auxiliary triple for any
two parity fields. It can recover the original rows through
a quartic double cover, with only the shared parity kernel
that vanished from the auxiliary support to be supplied.
That record is cheap on generic full-scale pairs and in
near-maximal fixed-target fibres. The original common
absolute residue gives positive numerical height
for the auxiliary even when the complete targets differ.
In the same-field case the two square-phase auxiliaries
also extend to moving and even targets with tall output
targets at this original absolute lock. Their actual
surviving prime orders identify when the strict
powerful-triple estimates are applicable.
[One-auxiliary determinant recovery (FP1–FP4)](ABC_FIELD_FREE_CROSS_PRODUCT_TRANSFERS.md)
· [Moving common-field transfer (CP6)](ABC_CROSS_PHASE_AND_AUXILIARY_TRIPLES.md#the-common-field-transfer-for-moving-and-even-targets)

## The unrestricted small-budget obligation

In the simultaneous packet sector, a sufficient missing
estimate is a fixed loss \(\sigma>0\) from
\(n^{2(1/6-d)}\) in a joint affine cell count with all three
outside radical caps. One can instead seek a loss from
\(n^{s_i+s_j}\) in the number of occupied pairs of absolute
residues for some projected role in each fixed packet tuple.
These are uniform bounds in the packet and residue data,
with the tolerance chosen after the loss and the height
taken to infinity at each fixed budget. A universal
argument must also handle the balanced profiles where this
near-minimum complete-dictionary reduction or the required
packet extraction has not been obtained.

In the simultaneous packet sector, bounding genuine shared-value
boxes below \(n^{md+(3m-1)t}\) by a fixed power for one fixed
\(m\) would now suffice. The known cycle and secant estimates
force large fresh prime supports. To turn them into an
exclusion requires an arithmetic upper bound on those
supports, or an independent upper bound for the arrays with
their full original caps. Generic private supports and full
power-class rank constrain a potential compression in this
sector. When \(w>6d\), ordinary synchronized differencing
cannot supply a second family at the same relative budget
on the forced Latin cells.

The remaining target is an absolute \(\delta>0\) and
\(N_\lambda(X)\ll_{\lambda,\varepsilon}
X^{(1/3-\delta)\lambda+\varepsilon}\) for all sufficiently small
fixed \(\lambda>0\). The least-target-repair equivalence retains
this obligation in the primitive quotient count. The graph
bounds give an exponent saving where a dense affordable affine
graph carries a small off-grid defect. The transfer bounds give
an exponent saving where a sufficiently imbalanced output has
bounded recorded radical and a numerically short remaining
cofactor. A proof at arbitrary profiles needs an estimate
forcing one of the requisite constraints on near-full input
families, or a different inequality using all three radical
caps. Families with broad positive valuation spectra on the
original letters, broad affine power-class distributions, and
large fresh radical in imbalanced outputs remain possible
within these necessary conditions. In the repaired formulation,
a saturated count must spread its representations over at least
\(Y^{(1/2-o(1))/M}\) occurring target radicals. Closing this
route requires a joint bound on how many supports have each
representation multiplicity, or the corresponding repeated
modular popularity estimate at block exponents \(\beta>1-u\)
when a short real window is imposed in the balanced normalization.
Equivalently it suffices to bound the mean capped signed test of
(CPH1) by \(N^{1-u-\sigma+\varepsilon}\) on a fixed finite mesh
covering the joint multiplicity range, with one absolute
\(\sigma>0\) and all sufficiently small fixed budgets.
Use the ramp with \(H=1/2\) at the occurrence endpoint \(u=0\). The corresponding valuation-class
route requires an entropy upper bound on simultaneous classes
with large effective least common multiples. The existing
powerful-number and Pell recovery estimates supply fibre
bounds after the relevant arithmetic records have been fixed.
They do not give either required distribution estimate.


In particular, an explicit sufficient extraction for closure is the following.
Fix \(k\) and one polynomial monomial support of invisible rank \(r<k\).
Choose absolute \(K,\gamma>0\) with
\[
                    k/K+\gamma<(k-r)/6 .
\]
In every putative balanced saturating sequence, find
\(n^{k/3-o(1)}\) ordered \(k\)-tuples for which these
least-repaired heights obey a nondegenerate relation
with that support, whose projective coefficient vector
runs through at most \(n^{\gamma+o(1)}\) records.
Equation (HP5) would contradict the tuple count,
uniformly in the original valuation spectra. The coefficient
entropy bound shows exactly how much compression this
route would have to extract from the original triples.
Coincidences involving only equal repaired targets,
and original additive relations with unbudgeted
repair ratios, do not supply the required tuples.
[Primitive lifting of repaired-height incidences](ABC_REPAIRED_HEIGHT_POLYNOMIAL_RELATIONS.md#transfer-back-to-primitive-abc-rows)
