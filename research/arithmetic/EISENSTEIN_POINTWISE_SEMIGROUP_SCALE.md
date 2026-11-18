# Pointwise semigroup-scale counts for balanced Eisenstein sums

Let \(A\) consist of \(1\) and the positive integers whose prime factors are
all \(1\bmod 3\). For \(n\equiv2\bmod6\), let \(r_H(n)\) count ordered pairs
\(a,b\in A\) with \(a+b=n\) and \(|a-b|\le H\). For \(n\equiv4\bmod6\),
let it count ordered pairs of actual summands \((a,3b)\), \(a,b\in A\), with
\(a+3b=n\) and \(|a-3b|\le H\). One may instead prescribe the orientation
\((3b,a)\). Define the finite product
\[
 J(n)=\prod_{\substack{p\mid n,\ p\ge5\\p\equiv2\bmod3}}
                  \frac{p-1}{p-2}.
\]

**Theorem EisSemigroupScale.** Fix \(2/3<\theta<1\) and put \(H=n^\theta\).
There exist positive constants \(c_\theta,C_\theta,N_\theta\) such that for
every integer \(n\ge N_\theta\) in either residue class,
\[
 c_\theta\,\frac{H J(n)}{\log n}
 \ \le\ r_H(n)\ \le\
 C_\theta\,\frac{H J(n)}{\log n}.
 \tag{1}
\]
In the mixed residue class each prescribed orientation satisfies (1).
The theorem determines the order of magnitude for individual targets. The
constants and onset are given existentially. An asymptotic with a specified
leading constant remains open here.

The separate [moving-window theorem EisLocationWindows](EISENSTEIN_MOVING_WINDOWS.md)
allows an arbitrary center \(x\), with matching count \(HJ(n)/\log n\)
whenever \(H\le M=\min(x,n-x)\) and
\(H^3\ge n^{1+\varepsilon}M\), for fixed \(\varepsilon>0\).
It retains the central threshold \(\theta>2/3\), and reaches designated
small-summand intervals of size \(n^\theta\) for every \(\theta>1/2\).

The lower bound uses the weighted shell of
[EisBoundedFactors](EISENSTEIN_BOUNDED_PRIME_FACTORS.md).
Requiring both components to have only large split prime factors gives a
sieve product of order \((\log n)^{-2}\). We allow three squarefree
cofactors made of small split primes. Two occur separately in the two
summands. The third is a common divisor of the summands and of the target.
Their combined contribution supplies the missing factor of \(\log n\)
uniformly, including when the target itself has many small split factors.

## The coordinate count and the cofactors

Work in \(R=\mathbf Z[\omega]\), \(\omega^2+\omega+1=0\), with norm
\(N(x)=x\bar x\). Call primes \(1\bmod3\) split and primes \(2\bmod3\)
inert. Put \(\delta=1-\theta\), \(h=H/n=n^{-\delta}\), and
\(b_h(t)=\max(1-|t|/h,0)\). In a shell \(N(x)+N(y)=m\) write
\[
 A_1=N(x),\qquad A_2=N(y),\qquad t=(A_1-A_2)/m .
\]
Count signed, ordered coordinates \((x,y)\) with weight \(b_h(t)\).
Its support excludes zero components for all sufficiently large \(n\).
Require both component norms to be odd. For mixed targets require also
\(9\nmid A_1,A_2\), meaning that neither is divisible by \(9\).
For each even \(m\) prime to \(3\) the odd shell has size
\[
 R_o(m)=12\{\sigma(m)-\sigma(m/4)\},
\]
where the second term is zero when \(4\nmid m\). For later normalization set
\[
 X_m=
 \begin{cases}
    (h/2)R_o(m),&m\equiv2\bmod6,\\
    (h/3)R_o(m),&m\equiv4\bmod6.
 \end{cases}
\]
The constant coefficient of the triangle is \(h/2\). The local mean of the
mixed condition at \(3\) is \(2/3\).

Take \(D=n^\xi,\ z=n^\eta,\ w=n^\lambda\) with positive exponents
chosen below and \(w\le z\). Consider pairwise coprime, squarefree positive
integers \(d,u,v\) made from split primes at most \(w\), subject to
\[
 d\mid n,\qquad (uv,n)=1,\qquad duv\le n^\kappa .
 \tag{2}
\]
For each triple work on the shell \(m=n/d\). Since \(d\equiv1\bmod6\),
\(m\) has the same residue modulo \(6\) as \(n\). Retain the original
relative width \(h\). After multiplying both component norms by \(d\), the
absolute width is \(dmh=H\).

Let \(W_{u,v}\) require \(v_p(A_1)=1\) at every \(p\mid u\), and
\(v_p(A_2)=1\) at every \(p\mid v\), where \(v_p\) denotes the prime
valuation. The opposite component is a \(p\)-unit because \(p\nmid m\).
At every other prime \(5\le p\le z\) require \(p\nmid A_1 A_2\).
This includes \(p\mid d\). Divide the coordinate weight by
\(2^{\omega(uv)}\), where \(\omega\) counts distinct prime factors.
This division compensates for the two choices of Eisenstein prime above
each prescribed split prime of valuation one.

## Exact divisibility at a split prime

At a split prime \(p\nmid m\), the two components of the Eisenstein
subalgebra modulo \(p^2\) occur on the diagonal of the split quaternion
matrix. The two coordinates of the other Eisenstein summand are the
off-diagonal entries, up to signs. The individual norm of either
summand is the product of its two entries, up to a unit. The finite
probability space is \(\mathrm{PGL}_2(\mathbf Z/p^2\mathbf Z)\) with the
determinant squareclass fixed by \(m\).

For reduction modulo \(p\), a fixed entry vanishes with probability
\(1/(p+1)\), and both diagonal entries vanish with probability
\(1/[p(p+1)]\). Indeed for determinant \(m\bmod p\) the first row is
uniform on nonzero rows and each row has \(p\) completions. There are
\(p-1\) matrices having both diagonal entries zero among the
\(p(p^2-1)\) matrices with that determinant. The same counts hold for
every nonzero determinant and for the off-diagonal entries. The only
compatible pairs among the four single-entry zero conditions are those
two pairs. Hence the probability of \(p\nmid A_1A_2\) is
\[
 B_p=1-\frac4{p+1}+\frac2{p(p+1)}
      =\frac{(p-1)(p-2)}{p(p+1)} .
\]
For \(v_p(A_1)=1\), choose which diagonal entry is zero modulo \(p\)
and require the other to be nonzero. The chosen entry must lift to
valuation one. Its conditional probability of such a lift is \((p-1)/p\).
Every matrix lift of an invertible residue is invertible, and unit
squareclasses modulo \(p^2\) for odd \(p\) are already detected modulo
\(p\). Scalar projectivization preserves these events. Thus the local
probability and its contribution relative to the rough sieve are
\[
 t_p=\frac{2(p-1)^2}{p^2(p+1)},\qquad
 a_p=\frac{t_p}{2B_p}=\frac{p-1}{p(p-2)} .
 \tag{3}
\]
The second component gives the same probability.

## The joint level estimate

**Lemma (joint prime-square level).** For squarefree \(\ell\) prime to \(6uv\), temporarily require
\(p\mid A_1A_2\) at all \(p\mid\ell\), keeping the prescribed factors,
oddness and the mixed condition. Before sifting, its weighted count on
the shell \(m\), including \(2^{-\omega(uv)}\), is
\[
 X_m\prod_{p\mid uv}(t_p/2)\,g_m(\ell)
 +O_\epsilon\!\left((uv)^3\ell^{3/2}n^{1/2+\epsilon}h^{-1/2}\right).
 \tag{4}
\]
Here \(g_m\) is multiplicative on squarefree integers and is the local
probability of \(p\mid A_1A_2\). We detail the extension of the
[composite level estimate](EISENSTEIN_BOUNDED_PRIME_FACTORS.md#the-composite-level-estimate)
needed for (4).

First symmetrize the entire \(W_{u,v}\) to
\((W_{u,v}+W_{v,u})/2\). The order units preserve the two component norms
or exchange them simultaneously. Thus this test, times all the other
symmetric local conditions, is left order-unit invariant. Multiplication
by the order unit \(j\) exchanges the two component norms. It preserves
the even triangle, the other divisibility tests and the shell. Consequently
the raw weighted count of this average equals the original ordered count
for each fixed \(u,v\). In computing the main term, apply the product
measure to each of the two products and take their average. Both product
means are equal. Separate swaps at different primes would describe a
different condition.

At all primes of the finite level coprime to the current target, use the
coefficient-valued all-ideal operator of
[EisPower](EISENSTEIN_POWER_SHRINKING_WINDOWS.md#uniform-archimedean-weight-and-finite-level).
It has norm at most \(\tau(m)\sqrt m\) at positive archimedean weight and
on the scalar cuspidal complement. On the remaining scalar norm characters
the main term is the product average conditioned on the target's
determinant squareclasses. These facts hold at any of the present product
levels by the class number one identification and the argument of
EisBoundedFactors (2). Adding the principal reduction level modulo
\(p^2\) at \(p\mid uv\) leaves the Hecke index coprime to the added
level. The unit norm characters at these odd primes still factor through
the determinant squareclass modulo \(p\), since scalar norms contain all
unit squares. The order of \(\mathrm{PGL}_2(\mathbf Z/p^2\mathbf Z)\) is
\(p^3|\mathrm{PGL}_2(\mathbf F_p)|<p^6\). Thus the added factors in
the finite Hilbert norm of a test bounded by one cost at most
\((uv)^3\). The squarefree coprime part of the old level costs at most
\(\ell^{3/2}\). Evaluation has norm at most one. Summing the positive
degree bounds against the absolutely convergent Legendre expansion of the
triangle costs \(O(h^{-1/2})\).

For completeness, the local probabilities at primes dividing \(\ell\)
are as follows. Write \(e=v_p(m)\). At split \(p\) with \(e=0\) the good
probability \(1-g_m(p)\) is \(B_p\). At split \(p\) with \(e\ge1\),
\[
 1-g_m(p)=\frac{(p-1)^2p^{e-1}}{(p+1)\sigma(p^e)} .
 \tag{5}
\]
At inert \(p\) the bad probability is
\[
 g_m(p)=
 \begin{cases}
  2/[p(p-1)], & e=0,\\
  0, & e=1,\\
  \sigma(p^{e-2})/\sigma(p^e), & e\ge2 .
 \end{cases}
 \tag{6}
\]
For inert \(p\mid m\), a bad point scales in both coordinates by \(p\)
from the shell at \(m/p^2\), or the bad set is empty. For split \(p\mid m\),
the [singular recurrence](EISENSTEIN_BOUNDED_PRIME_FACTORS.md#singular-primes)
stratifies by the minimum of the two component norm valuations at \(p\).
Its inversion replaces the bad count by angular counts at \(m/p^j\),
\(1\le j\le e\), unrestricted at \(p\), together with a bounded
modulo-\(p\) test at \(m/p^e\). Its coefficient absolute sum is at most
\(8e+2\). These transformations multiply both norms by the same
\(p\)-power and leave \(t\) unchanged. In particular they preserve the
two exact valuations at every other prime where such a prescription is
made. The extra test remaining at a now coprime split prime depends on
divisibility of the individual Eisenstein coordinates by primes above
that prime. Allocations of factors at distinct split primes in an
Eisenstein coordinate are independent by unique factorization. Thus the
recurrence applies to each of the two summands of the symmetrized test.
Every \(p\mid uv\) remains coprime to every reduced target. The joint
coefficient loss is at most \(\tau(m)^4\), as in EisBoundedFactors (6).
Finally, subtracting the even shell at \(m/4\) preserves the \(p^2\)
tests, since multiplication by \(2\) preserves their valuations. These
steps give (4), using \(\tau(m)^5\ll_\epsilon n^\epsilon\). They do not
require a bound for Hecke indices divisible by a level prime.

## Sifting outside the prescribed cofactors

Define
\[
 V_m(z)=\prod_{5\le p\le z}(1-g_m(p)).
\]
Uniformly in \(m\), \(g_m(p)=4/p+O(p^{-2})\) at split primes and
\(g_m(p)=O(p^{-2})\) at inert primes, with \(\sup_{m,p}g_m(p)<1\).
Use the elementary [block lower sieve](EISENSTEIN_BOUNDED_PRIME_FACTORS.md#an-elementary-lower-sieve)
outside \(uv\). The proof and its level are recalled here because the
error must be summed over the cofactors.

Partition the primes \(5\le p\le z,\ p\nmid uv\) into blocks
\(z^{2^{-j-1}}<p\le z^{2^{-j}}\), \(j\ge0\). In each block
\(\sum g_m(p)\le C\) and \(v_j=\prod(1-g_m(p))\ge c>0\) with absolute
constants, by prime harmonic Mertens. Omitting further primes preserves
both bounds. Let \(U_j,L_j\) be the usual inclusion-exclusion polynomials
for avoiding the bad primes in the block, cut at orders \(k_j\) and
\(k_j-1\), where \(k_j=2\lceil(K+j)/2\rceil\) and \(K\) is a fixed
large integer. Bonferroni gives \(L_j\le\mathbf1_{\rm good}\le U_j\)
and \(U_j\ge0\). Telescoping the product gives the pointwise lower bound
\[
 \mathbf1_{\rm all\ good}\ \ge\
 \prod_j U_j-\sum_j(U_j-L_j)\prod_{i\ne j}U_i .
\]
Under independent bad events of probabilities \(g_m(p)\),
\(\mathbf E U_j\ge v_j\) and
\(\mathbf E(U_j-L_j)\le C^{k_j}/k_j!\).
Choose \(K\) so that \(\sum_j C^{k_j}/(c k_j!)\le1/2\). Then the
expectation of the lower polynomial is at least
\(\tfrac12\prod_j v_j\).

Set \(S=\sum_{j\ge0}k_j2^{-j}<\infty\) and choose \(\eta=\xi/S\).
Every expanded divisor satisfies \(\ell\le z^S=D\). Because the blocks
are disjoint, each coefficient in the lower polynomial has absolute
value at most one plus the number of nonempty blocks, hence
\(O(\log\log(3n))\). Applying (4) yields, for each triple (2), a lower
main term of at least
\[
 \tfrac12 X_m V_m(z)\prod_{p\mid uv}a_p ,
 \tag{7}
\]
with absolute remainder
\[
 \ll_\epsilon (uv)^3 D^{5/2} n^{1/2+\epsilon} h^{-1/2}\log\log(3n).
\]
The factors \(a_p\) in (7) replace precisely the omitted \(B_p\), since
\(p\mid uv\) implies \(p\nmid m\). There are at most
\(\sum_{j\le n^\kappa}\tau_3(j)\ll_\epsilon n^{\kappa+\epsilon}\)
ordered triples, by counting all decompositions of \(j=duv\). Choose
\(\kappa,\xi>0\) such that
\[
 4\kappa+\tfrac52\xi<\tfrac12(1-3\delta).
 \tag{8}
\]
This is possible since \(\delta<1/3\). The sum of the absolute
remainders is at most
\(n^{1/2+4\kappa+(5/2)\xi+O(\epsilon)}h^{-1/2}\log\log(3n)\).
By (8), with sufficiently small \(\epsilon>0\), this is
\(o(H/\log n)\).

## Uniform contribution from the target divisors

The normalization of the shell is
\[
 \frac{R_o(m)}m=18\prod_{\substack{p^e\parallel m\\p\ge5}}
                              \frac{\sigma(p^e)}{p^e}.
 \tag{9}
\]
Indeed \(m\) is even and prime to \(3\), and
\(\sigma(2^e)-\sigma(2^{e-2})=3\cdot2^{e-1}\) for \(e\ge1\), with the
missing term at \(e=1\) taken as zero.
Use the target-independent product
\[
 V_*(z)=
 \prod_{\substack{5\le p\le z\\p\ {\rm split}}}B_p
 \prod_{\substack{5\le p\le z\\p\ {\rm inert}}}b_p,\qquad
 b_p=1-\frac2{p(p-1)} .
\]
It satisfies \(V_*(z)\asymp(\log z)^{-2}\). This follows by taking
logarithms, using \(B_p=1-4/p+O(p^{-2})\) and prime harmonic Mertens in
the two reduced progressions modulo \(3\).

At a split \(p\le z\) with \(p^e\parallel m,\ e\ge1\), (5) gives
\[
 \frac{\sigma(p^e)}{p^e}(1-g_m(p))
 =\frac{(p-1)^2}{p(p+1)}=B_p E_p,\qquad E_p=\frac{p-1}{p-2}.
\]
At inert \(p\le z\) dividing \(m\), (6) makes the same product
\(1+1/p\), since
\(\sigma(p^e)-\sigma(p^{e-2})=p^e+p^{e-1}\) for \(e\ge2\).
Its ratio to \(b_p\) is exactly \((p-1)/(p-2)\), also when \(e=1\).
At inert \(p>z\) dividing \(m\), the factor in (9) is at least
\(1+1/p=b_p(p-1)/(p-2)\). The product of these \(b_p\) is bounded
below uniformly by the positive product over all inert \(p\ge5\).
At split \(p>z\), the factor in (9) is at least one. All inert target
valuations are shared by \(m\) and \(n\) because \(d\) is split.
Since every \(E_p\ge1\), discarding the split enhancement for \(w<p\le z\)
therefore gives
\[
 X_m V_m(z)\gg
 \frac{H}{d}\,J(n)V_*(z)
 \prod_{\substack{p\le w,\ p\ {\rm split}\\p\mid m}}E_p .
 \tag{10}
\]

Temporarily drop \(duv\le n^\kappa\). The sum of the remaining positive
cofactor weights in (7) and (10) is an Euler product. At split \(p\le w\)
prime to \(n\), the possibilities are absence, membership in \(u\), or
membership in \(v\), with weights \(1,a_p,a_p\). At \(p^e\parallel n\),
\(e\ge1\), the possibilities are absence or membership in \(d\), with
weights \(E_p,E'_p/p\) where \(E'_p=1\) for \(e=1\), and \(E'_p=E_p\)
for \(e\ge2\). The factor \(1/p\) comes from \(h m=H/d\).
Their total satisfies
\[
 E_p+\frac{E'_p}{p}\ \ge\ E_p+\frac1p
   =\frac{p^2-2}{p(p-2)}=1+2a_p .
 \tag{11}
\]
It follows that the unrestricted Euler sum is at least
\(\prod_{p\le w,\ p\ {\rm split}}(1+2a_p)\asymp\log w\), since
\(a_p=1/p+O(p^{-2})\).

Normalize the local weights to independent probabilities. For \(p\nmid n\)
the probability of \(p\mid duv\) is \(2a_p/(1+2a_p)\ll1/p\). For
\(p\mid n\) it is \((E'_p/p)/(E_p+E'_p/p)\le1/p\). Consequently
\[
 \mathbf E\log(duv)\ll\sum_{p\le w}\frac{\log p}{p}\le C'\log w
\]
for an absolute constant \(C'\), by the prime sum bound and partial
summation. Choose \(0<\lambda<\min(\eta,\kappa/(3C'))\).
Markov's inequality now leaves at least half of the unrestricted Euler
sum on the triples satisfying (2). Inserting (10) and (11) into (7),
summing and using (8), gives total coordinate weight
\[
 \gg_\theta HJ(n)\frac{\log w}{(\log z)^2}
 \gg_\theta \frac{HJ(n)}{\log n}.
 \tag{12}
\]

## Removing the remaining inert primes and counting value pairs

For each \(d\) fix \(\rho_d\in R\) of norm \(d\); multiply one
Eisenstein factor above each split prime of \(d\). Map the residual
coordinates to the shell at \(n\) by
\((x,y)\mapsto(\rho_d x,\rho_d y)\). This is injective for fixed \(d\)
and preserves the angle, oddness, mixed condition, and divisibility by all
inert primes greater than \(z\).

The two original component norms recover the entire triple. Among split
primes at most \(w\), those of valuation exactly one in both original
norms are the divisors of \(d\). Those occurring with valuation one
only in the first, or only in the second, are the divisors of \(u\),
or of \(v\), respectively. The exact tests and the sieve at every other
small prime enforce these alternatives, including that both residual
norms are units at \(p\mid d\). Hence the images of distinct triples
are disjoint. Here we use their original ordered conditions, which
have the same weighted counts as their symmetrizations in (4).

Since the divided weights are at most one, the total loss from all
remaining inert primes is bounded by the original-target
[central tail estimate EisCentralTail](EISENSTEIN_POWER_SHRINKING_WINDOWS.md#a-large-prime-tail-inside-the-central-window).
Its upper bound is
\[
 \ll R_o(n)M(n)\left\{\frac h z+\frac{1+\log n}{\sqrt n}\right\}
   =o(H/\log n).
\]
Here \(M(n)\le\max_{k\le n}\tau(k)=n^{o(1)}\),
\(R_o(n)/n=n^{o(1)}\), \(z=n^\eta\), and \(\delta<1/3\).
Thus (12) persists with no inert prime factors of either norm.
Eisenstein norms are \(0\) or \(1\bmod3\), so at \(n\equiv2\bmod6\)
both are \(1\bmod3\). At \(n\equiv4\bmod6\) exactly one is divisible
by \(3\), and the mixed condition makes its valuation exactly one.
Exchanging the components and the cofactors \(u,v\) preserves the
entire weighted family. Each mixed orientation receives half the weight.

Write \(A_1=uk,\ A_2=vj\). All split primes of \(k,j\) exceed
\(n^\eta\), so each integer has at most \(1/\eta\) split prime factors
counted with multiplicity. The number of signed Eisenstein representations
of \(s>0\) is \(6g(s)\), \(g=1*\chi_3\), with local factors
\(g(p^e)=e+1\) for split \(p\) and \(g(3^e)=1\).
The only possible nonsplit prime of \(k,j\) is the single factor \(3\)
in the mixed component. Since the primes in \(uv\) occur to exponent one,
the residual signed coordinate multiplicity divided by \(2^{\omega(uv)}\)
is exactly \(36g(k)g(j)\le36\cdot2^{2/\eta}\). The triangle is at most one,
and the original value pair uniquely determines its triple and residual
value pair. Dividing (12) by this bound proves the lower half of (1).

## Upper bound for the same singular factor

The first actual summand belongs to an interval of length \(H\) centered
at \(n/2\). At each inert prime \(p\ge5\) a qualifying pair excludes the
classes \(0,n\bmod p\) for that first summand, with one excluded class
when \(p\mid n\). Multiplication of one summand by \(3\) does not change
this restriction. We may count all integers in the interval subject to
these exclusions. Put \(\nu(p)=1\) for \(p\mid n\), and \(2\) otherwise,
and extend multiplicatively to squarefree products of such inert primes.
For a given squarefree \(\ell\), CRT and interval counting give
\(H\nu(\ell)/\ell+O(\nu(\ell))\) integers lying in the excluded classes
at every \(p\mid\ell\), including the possible interval endpoints.

Use the same upper block polynomials with local bad probabilities
\(\nu(p)/p\), cutoff \(z'=n^{c}\), and corresponding fixed block sum
\(S'\). The block constants here may be enlarged independently of those in the lower sieve. They are uniform because \(\nu(p)\le2\).
The probability error for an upper block satisfies
\(\mathbf EU_j-v_j\le\mathbf E(U_j-L_j)\). Hence
\(\prod_j\mathbf EU_j\ll\prod_jv_j\), by the same convergent sum of
relative factorial errors. The expansion of \(\prod_j U_j\) has
coefficients of modulus at most one. Choose \(c>0\) sufficiently small
that its divisor level \(D'=(z')^{S'}\le H^{1/2}\). The total counting
error is bounded by
\(\sum_{\ell\le D'}2^{\omega(\ell)}\ll D'\log(2D')=o(H/\log n)\).
The remaining product is
\[
 \prod_{\substack{5\le p\le z'\\p\ {\rm inert}}}(1-2/p)
 \prod_{\substack{5\le p\le z',\ p\mid n\\p\ {\rm inert}}}
                             \frac{p-1}{p-2}
 \ll \frac{J(n)}{\log n}
\]
by prime harmonic Mertens. This proves the upper half of (1).

**Corollary (fixed central widths).** Fix \(0<c\le1\) and take \(H=cn\).
For every sufficiently large admissible \(n\), in the ordinary case and
in each prescribed mixed orientation,
\[
 r_{cn}(n)\asymp_c \frac{cn\,J(n)}{\log n}.
\]
In particular this applies to the central window \(H=n/3\) and to the
unrestricted count \(H=n\), since any two positive summands of \(n\)
have difference less than \(n\).

For the lower bound use the triangle with fixed relative width
\(h_0=\min(c,1/2)\). Its support is contained in the requested window,
and its shell main term has size proportional to \(h_0 R_o(m)\).
In the level estimate (4) the angular cost \(h_0^{-1/2}\) is fixed.
Choose \(4\kappa+(5/2)\xi<1/2\) in (8), the inequality with
\(\delta=0\). Then the summed remainders are \(o(h_0n/\log n)\).
The cofactor argument (7)--(12) gives a lower bound of order
\(h_0nJ(n)/\log n\), and the large-inert-prime tail is
\(o(h_0n/\log n)\) by the same displayed tail estimate. Since
\(h_0\ge c/2\), this gives the stated lower bound. The upper block
sieve applies directly to the requested interval of length \(cn\),
with \(D'\le(cn)^{1/2}\), and gives the matching upper bound.

The separate [joint coset estimate](EISENSTEIN_BOUNDED_PRIME_FACTORS.md#all-target-valuations-and-joint-cofactor-levels)
now improves the intermediate distribution remainder used by this
cofactor method. It handles every target valuation with one divisor
factor and gives an aggregate error proportional to the product of the
cofactor and sieve budgets, up to fixed logarithmic factors. The proof
above remains valid with its original sufficient parameter choices.

## Dependencies and review

The mathematical dependencies are
\[
 \text{angular Hecke bound and central tail (EisPower)}
 \ \longrightarrow\ \text{composite-level sieve and singular densities
 (EisBoundedFactors)}
 \ \longrightarrow\ \text{(4), (7) and EisSemigroupScale}.
\]
The angular Hecke estimate uses the unitary quaternionic operator
normalization, the positive-weight cuspidality check, Jacquet--Langlands
transfer and the classical Ramanujan bound as documented in
[EisPower](EISENSTEIN_POWER_SHRINKING_WINDOWS.md#cuspidality-at-positive-archimedean-weight).
The finite product-level conditional projection and singular-prime
inversion are given in [EisBoundedFactors](EISENSTEIN_BOUNDED_PRIME_FACTORS.md).
Primary sources for the coefficient model and the prime sums are
[Greenberg--Voight, Section 2](https://jvoight.github.io/articles/orthohecke-240123.pdf)
and [Loughran--Matthiesen, Lemma 2.4](https://arxiv.org/html/1904.12845v2#S2).
The block sieve uses elementary inclusion-exclusion, prime harmonic
Mertens and interval counting.

For comparison, Dietmann--Elsholtz, *Sums of two numbers having only prime
factors congruent to one modulo four*, Banach Center Publications 118 (2019),
177--179, [unnumbered Theorem, p. 177](https://doi.org/10.4064/bc118-11)
([author PDF, p. 1](https://www.math.tugraz.at/~elsholtz/WWW/papers/greaves20.pdf)),
gives eventual unbalanced existence for the Gaussian restriction. That stated
theorem concerns modulus four. The result here concerns the Eisenstein pure
and mixed families and gives the pointwise order of the number of pairs in
central power windows. This statement comparison does not settle historical
priority in the wider binary norm literature.

A separate internal reviewer reconstructed the cofactor argument,
including the prime-square lift, the simultaneous symmetry, the
target-divisor compensation and the conversion to distinct pairs. The
result is internally reviewed. External mathematical review and priority
assessment remain outstanding.
