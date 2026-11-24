# Rough Goldbach counts in random semigroups and quadratic fields

The random basis threshold survives when every witness in a positive
proportion of targets is required to be a sum of two squarefree
semiprimes. This remains true above the one-third roughness threshold.
Thus the threshold and the divisor-overlap estimates can hold on
actual representation families whose prime-pair coefficient is zero.
The result concerns the witness families and imposes no assertion
about the other representations of those targets.

A separate discriminant average reproduces the full prime-factor
generating function, with relative error uniform in the selection
probability. It works over odd negative fundamental discriminants of
size \(X\ge N^{1/2+\varepsilon}\) using Burgess. At that scale the
uniform quaternionic existence proof supplies no compatible weighted
lower bound. The two main terms for the Liouville test can already be
computed on average over even targets, and they are nonzero.

## Roughness and coefficient extraction

Fix an even \(N\), and let \(\mathcal R_\eta(N)\) be the set of ordered
pairs \(a+b=N\) with \(a,b>1\), \(ab\) squarefree and every prime
factor \(>z=N^\eta\). Asking for \(ab\) squarefree is the same as
requiring squarefree summands and coprimality. Strict roughness matches
the parity test. The version with primes \(\ge z\) has the same degree
bounds because \(a,b<N\). Put
\[
 R_N=|\mathcal R_\eta(N)|,\qquad
 F_N(t)=\sum_{\mathcal R_\eta(N)}t^{\Omega(ab)},\qquad
 G_N=[t^2]F_N(t).
\]
Each pair has \(\Omega(ab)\) distinct prime generators. Independence
therefore gives the stated exact expectation. With selection
probability \(tq^{-\gamma}\), its expectation instead has the extra
factor \((ab)^{-\gamma}\). On any fixed balanced window
\(a,b\asymp N\), this factor is comparable to \(N^{-2\gamma}\)
in every degree. There is no diagonal pair in \(\mathcal R_\eta(N)\).

Write \(J_\eta=2(\lceil1/\eta\rceil-1)\). Then
\(2\le\Omega(ab)\le J_\eta\). Positivity gives
\[
 0\le F_N(t)-t^2G_N\le t^3(R_N-G_N),\qquad
 G_N\ge \frac{F_N(t)/t^2-tR_N}{1-t}\quad(0<t<1).
 \tag{1}
\]
At \(\eta=1/5\) the degrees run from 2 through at most 8.
The two-dimensional lower sieve in
[Random prime semigroups](RANDOM_SPLIT_PRIME_SUMSETS.md#rough-decompositions-with-controlled-supports)
produces \(\gg_\eta N/(\log N)^2\) balanced pairs for every fixed
\(0<\eta<1/4.267\), so in particular at \(1/5\) and \(9/40\).
Indeed the \(d\)-th divisibility count is the CRT main term plus
\(O(2^{\omega(d)})\) in a fixed odd progression. With \(D=N^\xi\),
\(4.267\eta<\xi<1\), the weighted remainder is
\(O(D(1+\log D)^{15})\). Removing square factors and common primes
costs \(O_\eta(N/z+\sqrt N)\), and replacing a weak cutoff by a
strict one costs \(O(N/z+1)\) if \(z\) is a prime.
The [DHR lower-bound input, equations (6.11)–(6.12)](sources/johnston-2510.10853v3.pdf)
is positive at this choice of \(\log D/\log z\).

At \(1/4\) the degree bound is 6, and at \(1/3\) it is 4.
The lower sieve just described supplies neither threshold.
At \(1/3\), even a dimension-two level approaching \(N\) would give
\(\log D/\log z\le 3\), below the positivity range of this lower sieve.
At \(1/2\) every surviving summand is prime, and existence at that
cutoff already requires producing a prime pair.

For later use, for every fixed \(\eta>0\) the dimension-two upper
sieve gives
\[
 R_N\ll_\eta \frac{N}{(\log N)^2}S(N),\qquad
 S(N)=\prod_{\substack{q\mid N\\q>2}}\left(1+\frac1{q-2}\right).
 \tag{2}
\]
To see this, stop at \(y=N^{\min(\eta/2,1/20)}\), use level \(N^{1/2}\)
in the dimension-two upper sieve and the same CRT remainder, and
drop the squarefree and coprime conditions. At an odd \(q\) the two
deleted roots coincide exactly when \(q\mid N\), giving the factor
\((1-1/q)/(1-2/q)\) in \(S(N)\). The product of the generic factors is
\(\asymp(\log y)^{-2}\) by Mertens. Moreover
\(\sum_{N\le x}S(N)^2\ll x\): expand \(S(N)^2\) as
\(\sum_{d\mid N}g(d)\) over odd squarefree \(d\), with
\(g(q)=(1+1/(q-2))^2-1\ge0\); then
\(\sum g(d)/d=\prod_{q>2}(1+g(q)/q)<\infty\).
[Diamond–Halberstam, higher-dimensional upper and lower sieves](https://www.sciencedirect.com/science/article/pii/0022314X88900467)

The direct lower bound from the rough sieve is
\(F_N(t)\ge L_Nt^{J_\eta}\) with \(L_N\le R_N\).
Inserting it in (1) gives a nonpositive right side whenever
\(J_\eta\ge3\). Even exact knowledge of \(R_N\) leaves
\(R_Nt^{J_\eta}\) compatible with these inputs. A different option
is interpolation. For \(J_\eta-1\) fixed distinct positive points,
Lagrange interpolation of \(F_N(t)/t^2\) expresses \(G_N\) as a
fixed linear combination of the point evaluations. An error
\(o(N/(\log N)^2)\) per evaluation transfers at the same scale, with
a constant depending on \(\eta\) and the nodes. The coefficients
alternate, so this requires more than lower bounds at every node.

## Random representations by semiprimes alone

Consider any family of unordered pairs from \(\mathcal R_\eta(N)\)
with \(a,b\asymp N\), fixed support size \(j\), and
\(m\gg N/(\log N)^2\) members. At this target select eligible primes
independently with probability \(q^{-\gamma}N^{-\beta}\),
where \(\gamma,\beta\ge0\). If
\[
 2\gamma+j\beta<1,\qquad \gamma+\beta/\eta<1,                 \tag{3}
\]
the probability of no representative in the family is
\(O(\exp(-c N^\delta))\) for every fixed
\[
 0<\delta<\min\{1-2\gamma-j\beta,\ \eta(1-\gamma)-\beta\}.
 \tag{4}
\]
Indeed the success expectation is \(\mu\asymp
mN^{-2\gamma-j\beta}\). For squarefree \(d\), let \(m(d)\) count the
pairs whose product is divisible by \(d\). Assign each prime of \(d\)
to one of the two summands. For each of at most \(2^j\) assignments,
the CRT fixes one residue class for the first summand modulo \(d\).
Hence \(m(d)\ll_j N/d+1\), and \(\sum_{d>1}m(d)\le 2^jm\).
A common divisor \(g>1\) increases the joint success probability
by \(g^\gamma N^{\beta\omega(g)}\). Put \(w=N^{-2\gamma-j\beta}\)
and \(\gamma'=\gamma+\beta/\eta\). Enlarging the pair sum to a sum
over common divisors and splitting at \(d=N\) gives
\[
 \Delta \ll_j w^2\left(N^2\sum_{z<d\le N}d^{\gamma'-2}
                   + m/w\right)
 \ll_{j,\gamma'} w^2N^2N^{-\eta(1-\gamma)+\beta}+\mu.
\]
For the first range use \(N^{\beta\omega(d)}\le d^{\beta/\eta}\).
For the second use \(m(d)=O_j(1)\) and the maximum possible gain
\(N^{2\gamma+j\beta}=1/w\). Janson's bound
\(\exp[-\mu^2/(2(\mu+\Delta))]\) now gives (4).
At \(1/5\) the lower sieve and pigeonholing make this uniform for
every even target under (3) with \(j=8\), without identifying the
degree that supplies the witnesses.
[Janson–Warnke, equation (2)](https://arxiv.org/abs/1406.1248)

The general [random \(k\)-fold theorem](RANDOM_PRIME_SEMIGROUP_BASIS_THRESHOLD.md)
uses a \(k-1\) dimensional sieve to produce \(\gg
N^{k-1}/(\log N)^k\) tuples and discards the square and gcd exceptions.
After assigning a shared divisor to the \(k\) components, it
eliminates the coordinate with the smallest assigned divisor.
This gives the bound \(D(d)\ll N^{k-1}d^{-1+1/k}\) and leaves
a power saving in the Janson overlap sum for \(\gamma<1-1/k\).
The first moment uses \(R_u\asymp N^k\) for the full product of
the odd components. Neither step tracks prime-factor counts; the support size is
used only as a bounded constant. At and above the threshold, the almost-sure Euler-product
count, the limiting convolution volume and conditional Harris
bounds use all semigroup elements, also without factor strata.
The [fixed-multiplier refinement](RANDOM_PRIME_SEMIGROUP_LOGARITHMIC_TRANSITIONS.md)
has pole order \(\lambda\) for probabilities \(\lambda q^{-\gamma}\)
and uses conditional cofactors to recover logarithmic mass.
In its proper-overlap estimate \(\lambda^{-\omega(d)}\) is treated
as a fixed constant using the rough-support bound. Its estimates
for fixed \(\lambda\) are therefore a distinct use of the
multiplier from coefficient extraction as \(t\) varies with \(N\).

**Semiprime witness theorem.** For each fixed \(0<\eta<1/2\) and
all sufficiently large \(x\), \(\gg_\eta x\) even \(N\in[x,2x]\)
have \(\gg_\eta N/(\log N)^2\) balanced members of
\(\mathcal R_\eta(N)\) with \(\Omega(a)=\Omega(b)=2\).
For \(2\gamma+4\beta<1\), these pairs alone give
stretched-exponential success. Indeed one can construct them at
any fixed stricter cutoff \(\eta'\in[\tfrac13,\tfrac12)\) with
\(\eta'>\eta\), as (3) then follows from \(2\gamma+4\beta<1\).
For \(2\gamma+4\beta\ge1\), (2) and the union bound give success
probability \(O_\eta(S(N)/(\log N)^2)=o(1)\) uniformly for balanced
semiprime pairs. To bound \(S(N)\), use the Mertens product for
\(q\le\log N\). The remaining primes dividing \(N\) contribute
a bounded factor, since there are at most \(\log N/\log\log N\)
of them and each contributes \(1+O(1/\log N)\).

For the positive-proportion claim, choose \(\eta< u<v<1/2\).
There are \(\gg_\eta x/\log x\) squarefree semiprimes \(p q\) in
\([0.55x,0.75x]\) with \(p\in[x^u,x^v]\): count primes \(q\)
in \([0.55x/p,0.75x/p]\) by the PNT and sum \(1/p\) over
that range. The primes \(q\) are larger than \(p\), so this count
is injective, and both exceed \((2x)^\eta\) for large \(x\).
All sums fall in \([x,2x]\), and the semiprimes are odd.
Pairs sharing a prime cost \(O(x^2/(2x)^\eta)\) by summing the square
of the counts of multiples, even over integer moduli. The first
moment of the remaining ordered counts \(B_N\) is therefore
\(\gg_\eta x^2/(\log x)^2\). Equation (2) gives
\(\sum B_N^2\ll_\eta x^3/(\log x)^4\).
Remove targets with \(B_N<c_\eta x/(\log x)^2\) for a small
constant \(c_\eta\), then apply Cauchy–Schwarz to the remaining
first moment. They number \(\gg_\eta x\). Division by two gives
the unordered family used in (3).
[Mertens estimates used in the fixed-dimension sieve](sources/rosser_schoenfeld_1962.pdf)

The same positive-proportion theorem holds inside the pure
Eisenstein core. Require both \(p\) and \(q\) to be \(1\bmod3\) in
the semiprime construction. The fixed-modulus PNT gives the
same order of magnitude and puts every produced target in
\(N\equiv2\bmod6\). Each squarefree split semiprime has primitive
representations by the Eisenstein norm \(Q\), by factorization
in \(\mathbb Z[\omega]\). Thus all the stochastic witnesses in
this version lie on the actual quaternion shell and pass the
core congruence filters. Equations (2)–(4) then apply to the
same subfamilies of ordered norm values, counting each
unordered value pair once for Janson.

This is an arithmetic obstruction to reading prime pairs off the
threshold. On these subfamilies the weighted mean at \(\gamma=0\)
is exactly \(B_Nt^4\), with ordered counts. The same overlap input
and the same subcritical power-law exponent apply. The obstruction
does not concern stronger information such as the actual degree-two
part of the support hypergraph.

In fact only degree two has entirely disjoint supports automatically:
at fixed \(N\) a prime \(p\) determines its partner \(N-p\).
If \(P_0(t)\) is the probability that no restricted representation
survives, independence for prime pairs and Harris for the failures
of the other pairs give
\[
 (1-t^2)^{G_N/2}\!\prod_{\text{other unordered pairs}}
                 (1-t^{\Omega(ab)})
 \ \le P_0(t)\le (1-t^2)^{G_N/2}.
\]
Hence \(-\log P_0(t)=(G_N/2)(-\log(1-t^2))+O(t^3R_N)\)
for \(0<t\le1/2\). In the test window
\(t=N^{-\beta}\), \(1/3<\beta<1/2\), the error is \(o(1)\)
even using \(R_N\le N\). A uniform probability of success bounded
below there would force \(G_N\gg N^{2\beta}\).
The available input (3) with the worst-case degree 8 lies outside
this window. With all-prime supports supplied as input, it would
reach the window.

## The signed sums at the one-third cutoff

Fix \(1/3\le\eta<1/2\), and write \(L_N=\sum_{\mathcal R_\eta(N)}
\lambda(a)\) and \(K_N=\sum_{\mathcal R_\eta(N)}\lambda(a)\lambda(b)\).
By symmetry \(\sum\lambda(b)=L_N\).
Each summand has at most two factors, hence
\(4G_N=R_N-2L_N+K_N\). Let \(A_{ij}(N)\) count ordered pairs with
\((\Omega(a),\Omega(b))=(i,j)\), so \(A_{12}=A_{21}\).
With \(C_\eta=\log((1-\eta)/\eta)\), one has the unconditional
dyadic averages for the four cells \(i,j\in\{1,2\}\) at every
fixed \(0<\eta<1/2\)
\[
 \sum_{\substack{x<N\le2x\\2\mid N}} A_{ij}(N)
 \ \sim\ \frac{3x^2}{2(\log x)^2}c_i c_j,\qquad
 (c_1,c_2)=(1,C_\eta).                                       \tag{5}
\]
At \(1/3\le\eta<1/2\) these are all the cells, so
\[
 \frac{\sum L_N}{\sum R_N}\longrightarrow
       \frac{C_\eta-1}{C_\eta+1},\qquad
 \frac{\sum K_N}{\sum R_N}\longrightarrow
       \left(\frac{C_\eta-1}{C_\eta+1}\right)^2.
 \tag{6}
\]
At \(1/3\) these limits are about \(-0.1812\) and \(0.0328\).
Either uniform \(o(R_N)\) assertion would contradict (6).

For (5), first replace the moving cutoff by \(y=x^\eta\) or
\(y=(2x)^\eta\), squeezing each nonnegative cell separately.
For fixed \(s>0\), the number of permitted primes up to \(sx\)
is asymptotic to \(sx/\log x\). For distinct semiprimes it is
\(\sum_{y<p\le\sqrt{sx}}(\pi(sx/p)-\pi(p))\).
The total subtracted term is \(O(x/(\log x)^2)\).
The PNT in the first term, uniformly at arguments at least
\(\sqrt{sx}\), and partial summation over \(p\) give the integral
\(\int_\eta^{1/2}du/[u(1-u)]=C_\eta\), for either cutoff. Thus the
scaled counting measures on \([0,2]\) converge weakly to \(c_i\,ds\)
after normalization by \(x/\log x\). Every permitted summand is
odd. The area \(1<s+t\le2\) in the positive quadrant is \(3/2\),
and its boundary has zero limiting product measure.
Removing common primes costs \(O(x^2/x^\eta)\), as before.
Equation (6) follows by weighting the four cells.

There is a one-sided criterion that avoids a signed two-point estimate.
Let \(M_N\) count the pairs with a semiprime first summand. Then
\[
 G_N=R_N-2M_N+A_{22}(N)\ge R_N-2M_N.                          \tag{7}
\]
The ratio \(\sum M_N/\sum R_N\) tends to
\(C_\eta/(1+C_\eta)<1/2\). A pointwise proof would need a sifted
lower bound and an upper bound for a marked cofactor sum with a
strictly smaller constant. Obtaining a pointwise margin is a separate estimate.
For all \(\eta\), another one-sided certificate is
\(G_N\ge 3R_N-F'_N(1)\). At \(1/3\) it is exactly (7).
At \(1/5\), and more generally for \(0<\eta<1/(1+e)\), this
certificate is negative on a positive proportion of targets in
every sufficiently large dyad. Indeed
\(3R_N-F'_N(1)\le A_{11}(N)-A_{22}(N)\) since degree three
has zero coefficient and all degrees above it have negative
coefficients. By (5), \(\sum(A_{22}-A_{11})\) has positive
main constant \(3(C_\eta^2-1)/2\) at scale
\(x^2/(\log x)^2\). The positive part of \(A_{22}-A_{11}\)
has the same second moment upper bound as \(R_N\), from (2).
Cauchy–Schwarz therefore gives the positive proportion.
This rules out an approach based only on a uniformly positive
first-factor-moment certificate at the starting cutoff.

## What the quadratic norm proofs control

[Pointwise Correlations of Quadratic Form Values and Splitting Semigroups](https://hexagonmath.org/2611.00672v1)
uses a shared class of multiplying primes in condition (C2) and
a rough model with cutoff \(\exp((\log X)^{2/5})\). The stated
comparison also imports weighted dilation graph bounds and the
all-conductor zero-free half-plane in §3. Those imports do not fix
the elementary failed hypothesis in this application. Multiplication
by a bank prime below \(N^\eta\) destroys the requested roughness
on a nonzero input. On the unrestricted squarefree support,
adjoining a new prime multiplies the weight by \(t\) and changes
the Liouville sign. Including the bank record in the weight, or
retaining the predicate that it is empty, still breaks (C2).
The physical graph states also require marked prime divisors in
each bank. The polynomially rough endpoints furnish no such state.
The pole exponents and local phase factors count a full splitting
semigroup. They are not estimates of its fixed-factor coefficients
at polynomial roughness. This diagnosis does not require deciding
the independent proof status of the two recent imports.
[Comparison hypotheses](../../results/pointwise-binary-splitting-semigroups/comparison.tex) ·
[Imported inputs](../../results/pointwise-binary-splitting-semigroups/inputs.tex)

[Exact Sumsets of Primitive Quadratic Norms](https://hexagonmath.org/2611.00673v1)
and the [five Heegner cases](PRIMITIVE_NORM_HEEGNER_CORE_SUMSETS.md)
count a quaternion shell with norm \(N\), controlling the
coordinate-primitive conditions at inert primes and using exact
dyadic and ramified filters. At an inert \(q\nmid N\), a zero
component norm modulo \(q\) forces both of its coordinates to
vanish modulo \(q\), and the bad density in a unit norm class is
\(2/(q(q-1))\). The large-prime argument divides a nonzero bad
component by \(q\) and sums the resulting \(q^{-2}\) tail.
Split primes can divide a primitive norm through one of two
isotropic lines. Their bad density is of order \(1/q\) instead.
The core sieve neither deletes those factors nor distinguishes a
single such factor from several. Split-prime descent multiplies
both summands, so it also drops the desired coprimality and preserves
no upper bound on the factor count.

In the [uniform prime-discriminant proof](PRIME_DISCRIMINANT_CORE_UNIFORM.md#a-shell-estimate-with-the-discriminant-free),
the principal class mass is \(w=6/(p-1)\), the full degree is
of order at least \(N\), and the complementary all-ideal operator
has norm at most \(\tau(N)\sqrt N\). The inert normalizer test has
an absolute error on that spectral scale per auxiliary prime.
The ramified test removes zero blocks. For the large inert primes
the binary lattice count is \(O(X/\sqrt p+\sqrt X)\) at norm bound
\(X\), so \(P=(N/p)^{1/4}\) balances the spectral and area errors.
After normalization the loss is bounded by
\(O(\max_{v\le N}\tau(v)p^{3/4}N^{-1/4})\) for \(N\ge p^2\).
This gives the effective \(p^{3+\varepsilon}\) bound from the
explicit divisor estimate. None of these steps places a
factor-sensitive weight on the surviving split norms. The
all-class argument is important: a split prime need not be a
principal norm.

The [Eisenstein composite-level theorem](EISENSTEIN_BOUNDED_PRIME_FACTORS.md)
does test split divisibility. Its strengthened joint estimate
(J1)–(J5) gives macroscopic summed level \(D<N^{1/2-\varepsilon}\)
for the stated supported tests and coefficients. The good-prime
sieve dimension for the two norm blocks is two. Even at \(\eta=1/5\)
that spectral level gives \(\log D/\log z<5/2\), short of the
dimension-two DHR lower-sieve range. It does give a lower bound
for smaller fixed powers as recorded there, and gives no lower
bound at \(1/3\). A one-prime union-bound argument cannot
replace that lower sieve.

There is a concrete way to mark the missing arithmetic. In the
pure Eisenstein branch \(N\equiv2\bmod6\), restrict to squarefree
coprime \(a=Q(z), b=Q(w)>1\) at \(N^{1/3}\) roughness. Each such
nonramified value has \(6\cdot2^{\Omega(a)}\) coordinate witnesses,
by unique factorization and the six units. To preserve the
unweighted value count, use \(2^{-\Omega(ab)}/36\) on the shell.
For a sign test one can instead use any strictly positive weight
and ask for a sign-determined witness. For each retained summand,
\[
 \lambda(a)=-1+2\!\sum_{N^{1/3}<p\le\sqrt a}\mathbf1_{p\mid a}.
\]
The sum counts exactly one mark on a semiprime and none on a prime.
The one-mark term tests \(p\mid Q(z)\), and the two-mark term tests
both \(p\mid Q(z)\) and \(r\mid Q(w)\), with the two square-root
cutoffs. For split marks with \(p=\pi\bar\pi\),
\(r=\varrho\bar\varrho\), orientation followed by
\(z=\pi z'\), \(w=\varrho w'\) reduces the joint test to
\[
 pQ(z')+rQ(w')=N.
\]
Squarefreeness makes the orientations disjoint for each marked
summand. Thus (7) is a sufficient **one-mark** estimate on the
rough shell. The exact two-point test additionally requires the
corresponding **two-mark** estimate. These are explicit cofactor
families with \(N^{1/3}<p,r\le\sqrt N\).
The joint moduli \(pr>N^{2/3}\) lie beyond the proved summed
level. Even summing the available \(N^{1/2+o(1)}\) spectral error
over all the individual marks up to \(\sqrt N\) gives no
logarithmic saving at the target scale. The norm-shell Hecke
operator is indexed by the total norm \(N\). Its spectral bound
does not estimate signed sums across those varying cofactor
tests. A new distribution result for the rough shell and an
averaged estimate over the marked forms would be needed to apply
(7) in this branch. The [one-mark comparison](EISENSTEIN_ROUGH_ONE_MARK.md)
establishes the inequality for almost all targets by the circle
method and certifies the finite range through 20 million. A
pointwise lower bound for the signed minor-arc integral remains open.

## Discriminants as a source of prime selections

Let \(\mathcal D(X)=\{d:X<d\le2X,\ d\equiv3\bmod4,\ d
\text{ squarefree}\}\). The fields of discriminant \(-d\) have
local signs \(\chi_{-d}(q)=(-d/q)\) at odd primes \(q\). For an
odd squarefree \(m\le N^2\) with at most \(J\) prime factors
and any prescribed signs \(s_q\in\{-1,1\}\) at \(q\mid m\), their
joint probability over \(d\in\mathcal D(X)\) is
\[
 \prod_{q\mid m}\frac{q}{2(q+1)}+O_J(E),\qquad
 E=\frac{N^{1/2}\sqrt{\log(2N)}}{\sqrt X}.                      \tag{8}
\]
For every fixed \(\varepsilon>0\) there is also a power-saving
error \(E=O_{J,\varepsilon}(N^{-c_\varepsilon})\) uniformly
through \(X\ge N^{1/2+\varepsilon}\).
The first estimate is useful when one wants elementary conductor
dependence. It becomes a power saving at \(X\ge N^{1+\varepsilon}\).

For the proof, on \((d,m)=1\) expand the sign-pattern indicator as
\(2^{-\omega(m)}\sum_{e\mid m}(\prod_{q\mid e}s_q)(-d/e)\);
outside that coprimality condition the indicator is zero.
Expand \(\mu^2(d)=\sum_{h^2\mid d}\mu(h)\) and write \(d=h^2 u\).
Only odd \(h\) coprime to \(m\) contribute, with \(u\equiv3\bmod4\).
The \(e=1\) term has main count
\[
 \frac X4\frac{\varphi(m)}m
 \prod_{\ell\nmid 2m}(1-\ell^{-2}),
\]
with error \(O_J(\sqrt X)\) by counting the coprimality condition
with inclusion–exclusion for each \(h\le\sqrt{2X}\) and completing
the absolutely convergent sum. Dividing by the same formula at
\(m=1\) gives the product \(q/(q+1)\).
For \(e>1\) the character modulo \(e\) is primitive.
Expand the divisibility exclusions at \(m/e\) by inclusion–exclusion.
Each term is an interval character sum after removing a coprime
divisor and then using the progression modulo 4. Truncating the
sum over \(h\) at \(H\) gives Pólya–Vinogradov error
\(O_J(H\sqrt m\log(2m)+X/H+\sqrt X)\).
Use \(H\asymp\sqrt X/(m^{1/4}\sqrt{\log(2m)})\) when this is at
least 1 and the trivial estimate otherwise. This proves (8).
For the stronger range, use Burgess with an integer \(r>2\) and
any fixed \(\nu>0\) on each interval. Summation over all \(h\)
converges after extracting \(X^{1-1/r}\), since
\(\sum_{h\ge1}h^{-2+2/r}<\infty\). The normalized nonprincipal
error is \(O_{J,r,\nu}(X^{-1/r}m^{(r+1)/(4r^2)+\nu})\).
Squarefree \(e\) permits arbitrary \(r\), so take \(r\) large
in terms of \(\varepsilon\), then \(\nu\) sufficiently small.
[Heath-Brown, Burgess bounds, equation (1)](https://arxiv.org/abs/1203.5219)

**Uniform thinning theorem.** Draw \(j\) discriminants independently
from \(\mathcal D(X)\). At each eligible prime \(q\), choose a set
\(U_q\subseteq\{-1,1\}^j\) of accepted words, reject ramification,
and retain \(q\) when its word lies in \(U_q\). The mean of the
weighted representation count, where each retained pair carries
weight \(\prod_{q\mid ab}(1+1/q)^j\), agrees with the expectation in
the independent model with probabilities \(|U_q|2^{-j}\) up to
relative \(O_\eta(jE)\) if \(jE\) is sufficiently small. This bound
includes zero expectation, when both sides vanish.
In particular taking \(U_q=U\), of size \(t2^j\), gives
\(F_N(t)(1+O_\eta(jE))\). The raw uncorrected count has relative
error \(O_\eta(jE+jN^{-\eta})\) when this is small.

Indeed every allowed assignment of words on the support of a pair
has probability \(\big(\prod_{q\mid ab}q/(2(q+1))\big)^j
(1+O_\eta(jE))\). Independence here is across discriminants, and
(8) controls the signs within each discriminant. The main
probability in a single field is bounded below by a constant
depending only on \(J_\eta\). The number of word assignments is
\(\prod_{q\mid ab}|U_q|\), so summing their nonnegative contributions
preserves the relative error without a factor \(t^{-1}\).
One can take \(j=O(\log N)\) and the singleton word to realize
polynomially small \(t\) with a power-saving relative error. The
same calculation works for varying dyadic probabilities at each
prime. Interpolation only needs fixed \(j\). This theorem controls
the comparison with the independent model uniformly. It does
not provide a lower bound for either side at small \(t\).

For prime discriminants there is a simpler fixed-window comparison:
quadratic reciprocity expresses each sign as \((-1/q)(p/q)\).
CRT and Siegel–Walfisz give the fair moments over primes
\(p\in(X,2X]\) in a fixed class modulo 8 when \(N\) is
polylogarithmic in \(X\). The modulus for a chosen support divides
\(8ab\). The existing [growing-window theorem](RANDOM_SPLIT_PRIME_SUMSETS.md#transfer-to-a-power-sized-window-of-discriminants)
works up to \(N\le X^{1/12}\) for the **gap law** at density \(1/2\)
by majorizing a nonnegative tenth moment over primes with an
integer average. It gives a probability upper bound, not the
two-sided pointwise-in-\(N\) mean comparison in (8).

Heath-Brown's quadratic large sieve bounds the second moment over
squarefree conductors of \(\sum_{e\le N^2}c_e\chi_{-d}(e)\) with
odd squarefree \(e\), at cost
\((X+N^2)(XN)^\epsilon\sum_e|c_e|^2\).
The length \(N^2\) comes from the full support products \(ab\).
Such an upper bound does not give a positive small-\(t\) lower
bound or a signed improvement in the marked shell estimates. A
shorter discriminant average would need cancellation exploiting
the relation \(e\mid a(N-a)\) rather than a generic conductor
mean square.
[Heath-Brown, A mean value estimate for real character sums, Theorem 1](https://www.impan.pl/en/publishing-house/journals-and-series/acta-arithmetica/all/72/3/108575/a-mean-value-estimate-for-real-character-sums)

There is a concrete range mismatch. The proven uniform
prime-discriminant core bound guarantees existence when
\(N\ge C_\epsilon p^{3+\epsilon}\), whereas (8) with Burgess uses
\(X\ge N^{1/2+\varepsilon}\) and a different conductor family.
The principal representation obstruction is stronger in the
elementary range. If \(d>4N\), the identity
\(4Q_d(x,y)=(2x+y)^2+dy^2\) implies that every principal norm
at most \(N\) is a square. The restricted principal shell is
then empty, even though the splitting-semigroup average can be
positive. Averaging arbitrary ideal norms can recover splitting,
but the principal quaternion shell estimates do not give a
signed lower bound for that average.

Finally, \(\lambda(n)=\chi_D(n)\) on the semigroup generated
solely by inert primes, while \(\chi_D(n)=1\) on the split
semigroup. The inert coset can encode parity at a fixed conductor:
in the all-inert family at \(\eta\ge1/3\), requiring both
character values \(-1\) counts exactly its prime pairs. At
\(\eta=1/5\), it allows \(\Omega=1,3\) in each component, so
the contamination begins in degree four. This exact filter
doesn't restore the missing roughness estimates. The fixed
prime-coset comparison again uses prime multiplication at the
subpower bank and its exact invariance hypothesis.

## Comparison with prime-producing sieves

- Chen's prime-plus-\(P_2\) result uses a weighted sieve and
  switching. Its rough factor bound in the original announcement
  is \(N^{1/10}\) on the complementary part. The Li–Liu preprint
  gives \(p+rq\) with \(r\le q^{0.9}\) and a positive lower
  count. This restricts the imbalance but does not force \(r=1\)
  or \(r>N^{1/3}\). The November Slava Chan preprint claims
  the tighter \(r\le q^{0.75}\); its full proof was outside
  the present audit. These refinements are distinct from the
  composite-only random threshold above.
  [Chen's 1966 announcement](https://www.sciengine.com/doi/pdf/8663c8f4bb65481c80325a4101d215ac) ·
  [Li–Liu, Theorem 1.1](https://arxiv.org/abs/2606.05224) ·
  [Chan preprint](https://aixiv.science/pdf/aixiv.261111.000001v1.0.pdf)
- Friedlander–Iwaniec's asymptotic sieve adds a bilinear axiom
  used for primes of the form \(a^2+b^4\).
  Heath-Brown's \(x^3+2y^3\) theorem exploits multiplicative
  structure in cubic ideals to verify bilinear estimates.
  Ford–Maynard give both positive prime-detection criteria and
  examples that check the adequacy of prescribed Type I and II
  ranges. The quaternionic norm of a product is multiplicative
  for the **total** norm, which does not supply those hypotheses
  for either individual block along the fixed additive constraint.
  The explicit marked forms above are where comparable additional
  arithmetic would have to enter.
  [Friedlander–Iwaniec, Asymptotic sieve for primes](https://arxiv.org/abs/math/9811186) ·
  [Heath-Brown, Primes represented by \(x^3+2y^3\)](https://ora.ox.ac.uk/objects/uuid%3A8b6404f0-0646-4bc9-8dd4-c2ca735ab3b9) ·
  [Ford–Maynard, On the theory of prime-producing sieves](https://arxiv.org/abs/2407.14368)
- Lichtman and Pascadi improve levels for specified factorable
  weights. The Goldbach application in Lichtman is an upper
  bound. Those results do not supply the marked rough
  quaternionic distribution. Evans's restricted almost-prime
  correlation theorem averages over difference shifts with
  one factor size restricted. Tao's logarithmically averaged
  Chowla theorem and Mangerel's reflected Liouville results
  concern different averaging or support. In particular,
  a reflected unrestricted sign pattern includes products
  of three or more primes, and (6) rules out importing a
  zero-centered estimate onto the rough family.
  [Lichtman](https://arxiv.org/abs/2309.08522) ·
  [Pascadi](https://arxiv.org/abs/2505.00653) ·
  [Evans](https://arxiv.org/abs/2102.12297) ·
  [Tao](https://arxiv.org/abs/1509.05422) ·
  [Mangerel](https://arxiv.org/abs/2412.17199)
- The recent *A Goldbach theorem for the Liouville function*
  states unconditional all-even existence on the unrestricted
  domain. In §§2–3 the descent scales both summands by a shared
  factor, and the reflected sign propagation builds on doubling
  both parts of a sum and transporting signs by 2 and 3.
  On the rough domain the doubled parts are automatically excluded
  once \(N^{1/3}>2\). The missing-pair assumption there cannot
  start the propagation in Lemma 3.1. The scaling step also
  creates a common divisor. Hence importing the two commuting
  symmetries does not establish rough sign patterns.
  [A Goldbach theorem for the Liouville function, §§2–5](https://ggeta.com/papers/liouville-goldbach-gpt6-astra.pdf)

The positive contributions are the weighted discriminant transfer,
the actual composite-only witness theorem, and the centering formula
(6). Their proofs use classical inputs with the repository's support
and norm-shell setup. They improve neither the all-target Goldbach
lower bound nor Chen-type representability. No independent review or
priority determination is claimed.

## Finite experiments

Run `python3 research/arithmetic/experiment_rough_goldbach.py`
from the repository root. The standard-library script uses exact
integer roughness tests, least prime factors, exact \(\Omega\)
and coprimality. It records ordered coefficients and the four
parity cells. The discriminant calculation enumerates squarefree
\(d\equiv3\bmod4\), counts the splitting event for each support,
and raises its exact frequency to the \(j\)-th power to average
over independent fields. Relative errors are rounded only for
display.

| \(N\) | \(\eta\) | Ordered coefficients starting with degree 2 |
| ---: | ---: | :--- |
| 1,000 | \(1/5\) | \(54,58,16\) |
| 1,000 | \(1/3\) | \(54,38,4\) |
| 1,000,000 | \(1/5\) | \(10804,25966,20712,5788,548,2\) |
| 1,000,000 | \(9/40\) | \(10802,22654,14650,2798,118\) |
| 1,000,000 | \(1/3\) | \(10796,11662,3216\) |

Trailing zeros are omitted. The \(N=10^6,\eta=1/5\) row gives
\(3R_N-F'_N(1)=-23136\), although \(G_N=10804\).
For the requested extraction expression \(F_N(t)/t^2-tR_N\)
the same row gives \(1396521/512>0\) at \(t=1/4\) and
\(-34995/16<0\) at \(t=1/2\), using exact finite counts.
At \(1/3\) the
certificate \(R_N-2M_N\) is positive for each tested even target
from 100 through 5000. The \(N=10^6,\eta=1/3\) row gives
\(L_N/R_N=-0.2952\), \(K_N/R_N=0.0915\) and \(R_N-2M_N=7580\).
The composite-only family at this target has 3216 ordered witnesses,
all invisible to the prime-pair coefficient.

At \(N=1000,\eta=1/5\) and \(X=65536\), there are 13270 sampled
discriminants in the dyad. The ramification-corrected expectations
for \(t=1/2,1/4,1/16\) have relative errors \(+0.441\%\),
\(+0.600\%\), and \(+0.680\%\). The errors vary in sign with
the dyad. These finite observations check the identities and
illustrate the ranges. They supply no pointwise asymptotic or
uniform sign claim.
