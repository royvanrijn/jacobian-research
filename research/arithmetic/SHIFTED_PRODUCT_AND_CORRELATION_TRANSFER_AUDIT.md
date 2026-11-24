# Shifted powers and Gaussian norm triple extraction

## Shifted powers with a summable coefficient family

Let \(k\ge2\) be an integer, \(h\) a fixed nonzero integer, and
\((w_j)\) a nonempty finite or countable family of positive integers
satisfying \(\sum_j w_j^{-1/k}<\infty\). Define
\[
 M=\bigcup_j \{w_j x^k+h:x\in\mathbb Z_{\ge1}\}\cap\mathbb Z_{\ge1}.
\]
Suppose \(R\subset\mathbb Z_{\ge1}\) satisfies
\[
 |(R\setminus M)\cap[1,X]|=o(X^{1/k}),\qquad
 \limsup_{X\to\infty}X^{-1/k}|R\cap M\cap[1,X]|>0.
\]
Then \(R\) and every finite modification of \(R\) have no factorization
\(A B\) with \(A,B\subset\mathbb Z_{\ge1}\) and \(|A|,|B|\ge2\).
The error hypothesis \(|(R\mathbin\triangle M)\cap[1,X]|=o(X^{1/k})\)
is sufficient because one coefficient family already supplies
\(\gg X^{1/k}\) values. In particular \(w_j=j^\ell\) covers
\(x^k y^\ell+h\) for \(\ell>k\), including coprime exponents and
\(\min(k,\ell)=2\). Equal exponents give a single coefficient
because \(x^k y^k=(xy)^k\). For a finite product of powers, combine
all factors having the minimum exponent into one variable. The
remaining coefficient series is a finite product of convergent
reciprocal-power series.

Hajdu and Sárközy posed the single-family \(h=1\) problem in
Conjecture 1.2 of paper III and the product problem in Conjecture 1
of paper II. Yip's Theorem 1.3 proves the fixed-shift single-family
theorem for \(k\ge3\). The supplied catalogue proves both all-order
theorems and gives a stronger shared hypothesis than symmetric
difference: the excess \(|(R\setminus M)\cap[1,X]|\) is
\(o(X^{1/k})\) and the retained count divided by \(X^{1/k}\) has
positive upper limit. The determinant proof below independently
checks the source mechanism and uses a one-dyad sum on each
fixed macroscopic product band. Finite changes preserve either
hypothesis, which is the total primitivity convention. The
catalogue also supplies Lean artifacts. The argument here is a
written proof check.
[Hajdu–Sárközy II](https://math.unideb.hu/sites/default/files/inline-files/poly2_revised.pdf)
· [Hajdu–Sárközy III](https://real.mtak.hu/102769/1/poly3_190407.pdf)
· [Yip, arXiv:2501.16620v3](https://arxiv.org/html/2501.16620v3)
· [Catalogue shifted-power proof](https://github.com/1337isnot1337/solved-conjectures/blob/758066055ef49c2c66cb1087918ae62ad75e7ac1/results/hajdu-sarkozy-shifted-power-perturbations/PROOF.md)
· [Catalogue summable-family proof](https://github.com/1337isnot1337/solved-conjectures/blob/758066055ef49c2c66cb1087918ae62ad75e7ac1/results/hajdu-sarkozy-shifted-products/PROOF.md)

For distinct \(a_1<a_2\) and \(b_1<b_2\), put \(s_{ij}=a_i b_j-h\).
Then
\[
 s_{11}s_{22}-s_{12}s_{21}=-h(a_2-a_1)(b_2-b_1).             \tag{1}
\]
If all products \(a_i b_j\) lie in \([Y,(1+\delta)Y]\), then
\((a_2-a_1)(b_2-b_1)\le\delta^2 a_1 b_1\le\delta^2(1+\delta)Y\).
If all \(s_{ij}=w x_{ij}^k\) for the same \(w\), the nonzero
left side of (1) has absolute value at least
\[
 k w^{2/k}(Y-|h|)^{2-2/k}.
\]
This follows by applying the mean value theorem to the two distinct
positive integer products \(x_{11}x_{22}\) and \(x_{12}x_{21}\).
Since \(w\ge1\) and \(2-2/k\ge1\), a fixed sufficiently small
\(\delta>0\) depending on \(k,h\) excludes every four-cycle in that
product range for all sufficiently large \(Y\), uniformly in \(w\).

A second use of (1) fixes \(a_1<a_2\), allows one coefficient \(w_1\)
on the first row and another \(w_2\) on the second row, and restricts
\(b_1,b_2\) to \([Z,(1+\delta')Z]\). The right side is at most
\(|h|(a_2-a_1)\delta' Z\). For large \(Z\), the left side has absolute
value at least
\[
 k (w_1w_2)^{1/k}(a_1 a_2 Z^2/4)^{1-1/k}.
\]
Indeed each product \(s_{11}s_{22}\) and \(s_{12}s_{21}\) is at least
\(a_1a_2 Z^2/4\), and they are distinct \(k\)-th powers multiplied
by the common coefficient \(w_1w_2\). Choose \(\delta'\) small
depending on \(a_1,a_2,k,h\). There is at most one simultaneous
neighbor per interval in the resulting geometric cover, hence
\(O(\log(2X))\) neighbors up to \(X\) for each pair of coefficients.

Set \(\rho=1/k\) and suppose \(R=AB\). If \(M_J\) is the union of
the first \(J\) coefficient families, then
\[
 |(M\setminus M_J)\cap[1,X]|
 \le (X+|h|)^\rho\sum_{j>J}w_j^{-\rho}.
\]
Fix two distinct \(a_1,a_2\in A\). Every \(b\in B\) up to \(X\) has
\(a_i b\in R\). If both products belong to \(M_J\), the preceding
neighbor bound gives \(O_J(\log(2X))\) choices. The other choices are
bounded by \(\sum_{i=1}^2 |(R\setminus M_J)\cap[1,a_i X]|\), since
each map \(b\mapsto a_i b\) is injective. Divide by \(X^\rho\), take
an upper limit, and let \(J\) increase. This gives \(B(X)=o(X^\rho)\).
Using two fixed elements of \(B\) gives \(A(X)=o(X^\rho)\).

Fix \(c\in(0,1)\). Cover \([cX,X]\) by a bounded number of the
multiplicatively thin intervals above. For each product orientation
we may take \(a\le\sqrt X\), and dyadically \(D\le a<2D\), \(D\le\sqrt
X\). At fixed \(w_j\) in one thin interval, (1) allows at most one
common neighbor to two distinct left vertices. If the two vertex
counts are \(P,Q\) and their edge count is \(E\), summing
\(\binom{\deg b}{2}\) over the \(Q\) right vertices gives
\(E^2/Q-E\le P(P-1)\), hence \(E\le P\sqrt Q+Q\). Take
\(P\le A(2D)\), \(Q\le B(X/D)\). The sum of the first terms is
\[
 \ll X^{\rho/2}\sum_{D\le\sqrt X,\ D=2^j} D^{\rho/2}
 \ll X^{3\rho/4}=o(X^\rho)
\]
by the upper bounds \(A(U),B(U)\ll U^\rho\). For the second terms
\(B(X/D)=o((X/D)^\rho)\) uniformly since \(X/D\ge\sqrt X\), and
\(\sum_D D^{-\rho}\ll1\). The other orientation has the same bound
with \(A,B\) exchanged. Summing the finitely many terms for \(M_J\)
therefore counts at most \(o_J(X^\rho)\) represented values of
\(M_J\) in the band. Take \(c=1/2\) and sum these bounds over the
dyadic bands below \(X\). Given \(\eta>0\), all sufficiently large
upper endpoints \(Y\) contribute at most \(\eta Y^\rho\), and the
sum of these powers is \(O_\rho(X^\rho)\). The bounded small range
has negligible count, giving \(|(AB)\cap M_J\cap[1,X]|=o_J(X^\rho)\).
The tail bound shows that the normalized upper limit for
\((AB)\cap M\) is at most \(\sum_{j>J}w_j^{-\rho}\). Let \(J\)
increase to contradict the positive retained upper limit. Finite
changes preserve both hypotheses.

For an additive rectangle \(u_{ij}=a_i+b_j\in[Y,2Y]\) of squares,
the analogous local exclusion has a different scale. The nonzero
integer \(\sqrt{u_{11}}+\sqrt{u_{22}}-\sqrt{u_{12}}-\sqrt{u_{21}}\)
is the integral of \(-1/(4t^{3/2})\) over the two side directions.
Thus \((a_2-a_1)(b_2-b_1)\ge4Y^{3/2}\). In particular a pair of
rows at gap \(d>0\) has at most \(d/(2\sqrt Y)\) common square
neighbors on this band. Each such neighbor gives a positive integer
difference \(r-s\le d/(2\sqrt Y)\), and \(d=(r-s)(r+s)\) determines
the pair once the difference is fixed. These two constraints hold
without a divisor-function loss.

The supplied logarithmic-loss theorem uses a different global
argument, which survives an independent check. For \(R=A+B\), let
\(e_\pm(X)\) count missing squares and added nonsquares, write
\(e=e_-+e_+\), and use \(K\) for the Landau–Ramanujan constant.
The hypothesis \(e(X)\le\kappa\sqrt{X/\log X}\) with fixed
\(\kappa<K\) leaves \(\gg X/\sqrt{\log X}\) different sums of two
kept squares whose crossed coordinate swap is still a pair of
squares. Missing squares and failed swaps remove in total at most
\(e(X)(\lfloor\sqrt X\rfloor+1)+e_+(X)^2\) values. In the remaining
swaps, a fixed ordered pair of distinct left coordinates with
gap \(d\) gives at most \(\tau(d)^2\) sums by factoring the two
equal square differences. For repeated left coordinates, use
the right gap, with a separate \(O(\sqrt X)\) diagonal. This
forces one factor size to be \(X^{1/2-o(1)}\). Square codegree
bounds and two fixed translates also give both factor sizes at
least \(X^{1/4-o(1)}\) and both \(o(\sqrt X)\).

In the isolation argument, each chosen square in \([X/2,X]\) has
an interval of radius \(\delta\sqrt X\) containing no other element
of \(R\). It has two coordinates isolated at the same distance
in \(A,B\). With one coordinate \(a\ge X/4\) and the smaller one
in \([Y,2Y)\), let \(m_Y,n_Y\) count the used vertices. Short
translations by the opposite set have disjoint supports and at
most one square per interval, giving
\[
 m_Y(B(\delta\sqrt X/2)-1)\le e(2X),\qquad
 n_Y(A(\delta\sqrt Y/2)-1)\le e(3Y).
\]
Their product gives \(m_Yn_Y\ll X^{1/2}Y^{1/8+o(1)}\) using
\(A(t)B(t)\ge t^{3/4-o(1)}\). The symmetric codegree bound gives
\(O(m_Y+n_Y+X^{o(1)}(m_Yn_Y)^{3/4})\) square edges. The linear
terms sum to \(o(\sqrt X)\) because the smaller-coordinate bins
are disjoint and \(B(\delta\sqrt X/2)\ge X^{1/8-o(1)}\). The
nonlinear terms sum geometrically to \(X^{15/32+o(1)}\).
Swapping the orientations covers all isolated squares and
contradicts their \(\gg\sqrt X\) count. The same proof reaches
finite modifications by leaving slack in the choice of \(\kappa\).
The square-swap argument concerns full sumset equality and
sparsity of its added values. The fixed-target square-full
representation graph supplies neither hypothesis.
[Catalogue perturbed-squares proof, §§1–3](https://github.com/1337isnot1337/solved-conjectures/blob/dab1ac8d1ecea8dc34c38932fcd2a4b3f56ce0ef/research/erdos-perturbed-squares/LOG-LOSS.md)
· [Sárközy–Szemerédi, Theorem II](https://real-j.mtak.hu/9363/1/MTA_MatematikaiLapok_1965.pdf#page=78)

## Determinants in the square-full kernel graph

For a coprime edge of the graph on squarefree kernels \((b,d)\), write
\(a^2b^3+c^2d^3=N\). On a cycle of length \(2r\), let \(U,V\) be
the two products that each take \(r\) left-root labels \(a\) and
\(r\) right-root labels \(c\), in alternating orientations. Multiplying
the congruences \(a^2b^3\equiv-c^2d^3\bmod N\) and cancelling the
kernel cubes gives \(N\mid U^2-V^2\). Every label is a unit
modulo \(N\) by coprimality of the representation. For any odd
prime power \(q\mid N\), one sign satisfies \(q\mid U\pm V\).
If \(a\le A\), \(c\le C\), then \(q>2(AC)^r\) forces \(U=V\).
A factor \(2^e\mid N\) gives the same implication with \(2^{e-1}\)
in place of \(q\), by the valuations of the sum and difference of
two odd integers. In particular, all four-cycles in the balanced
box are singular for sufficiently large prime \(N\), since there
\(A,C\ll N^{1/5}\). For six-cycles the available height is
\(O(N^{6/5})\).

The integer-gap argument for shifted powers used a determinant
whose small nonzero value was fixed by the two side gaps. The
kernel determinant can vanish while the kernels remain formally
distinct. For example, if positive integers \(p_1,p_2,q_1,q_2\)
parametrize the rank-one edge ratios \(p_i/q_j\), one recovers
coprime roots \(p_i/g_{ij},q_j/g_{ij}\) with \(g_{ij}=(p_i,q_j)\).
The alternating square-sum test for a singular four-cycle is
\(g_{11}^2+g_{22}^2=g_{12}^2+g_{21}^2\). The choices
\(p_1=4,p_2=56,q_1=7,q_2=8\) give gcd matrix
\(\left(\begin{smallmatrix}1&4\\7&8\end{smallmatrix}\right)\)
and \(1+64=16+49\). These are exact ratio and gcd constraints.
Constructing a cycle would also require positive squarefree kernel
cubes solving the edge equations. The singular and full-height
nonsingular cases of GAP-ARITH-ERDOS943-BALANCED-SIX-CYCLE-POINTWISE
remain open, as does GAP-ARITH-ERDOS943-POINTWISE-KERNEL-PAIR-SAVING.
[Cycle equations](ERDOS_943_SIX_CYCLES.md)
· [Kernel-pair bound](POWERFUL_SUM_GRAPHS.md#square-full-kernel-pairs-at-the-balanced-scale)

## Triple extraction from Gaussian norm pair moments

The supplied smooth-number argument produces the pattern \(010\),
with exactly one hit in the middle, and excludes any summand that
gives all its elements a neighbor at one prescribed distance.
Its proof uses stationarity and pair marginals. For a stationary
binary process with density \(0<d<1\) and pair probabilities
\(d^2\) at all distinct offsets, expansion at step size \(h>0\)
gives \(p_{001}=p_{010}\), since their common triple term cancels
from the difference. A length-\(H\) word starting \(00\) either
contains \(001\) at the first hit or is empty. Its count has mean
\(Hd\) and variance \(Hd(1-d)\), so
\[
 (1-d)^2\le (H-2)p_{010}+(1-d)/(Hd).
\]
Taking \(H>1/[d(1-d)]\) forces \(p_{010}>0\).
The proof remains valid with asymptotic first and pair marginals
on an asymptotically translation-invariant window. The catalogue
gets these for the smooth-number indicator in arbitrarily slowly
growing harmonic windows. Teräväinen's Theorem 1.4 applies to
the frozen cutoff \(X^\epsilon\) for \(0<\epsilon<1\). On bounded
progressions the relevant excluded products of primes
\(>X^\epsilon\) are units, have bounded numbers of factors, and
number \(O_\epsilon(X/\log X)\), so inclusion–exclusion gives
the required \(o(X)\) progression discrepancies. Dickman
asymptotics give the dyadic mean. Mertens and
\(X^{\epsilon-1}\omega(X)\to0\) bound the cost of replacing the
frozen cutoff by \(n^\epsilon\) on
\([X/\omega(X),X]\), \(1\ll\omega(X)\le\log(3X)\).
If the natural lower density of a pattern were zero, a subsequence
with density \(D_X\to0\) and an admissible \(\omega(X)\) chosen
with \(\omega(X)D_X\to0\) would contradict the harmonic lower
bound. The separate finite arbitrary-factor extension assumes
Schinzel's hypothesis H and does not enter this deduction.
[Catalogue fixed-gap note](https://github.com/1337isnot1337/solved-conjectures/blob/dab1ac8d1ecea8dc34c38932fcd2a4b3f56ce0ef/research/sarkozy-smooth/FIXED-GAP.md)
· [Teräväinen, Theorem 1.4 and §4](https://arxiv.org/pdf/1710.01195)
· [Győry–Hajdu–Sárközy, Part I](https://arxiv.org/pdf/2006.15307)

The fixed-gap result concerns forbidden patterns of sumsets with a
prescribed neighbor property. The translation-pattern goal
ARITH-ERDOS675-SUMS-TWO-SQUARES-TRANSLATION asks for a repeat of
each individual finite initial word. A positive density for an
individual short word leaves that recurrence obligation open.
[Exact translation goal](ERDOS_675_TRANSLATION_PATTERNS.md)

The smooth-number density \(d\) is fixed. For Gaussian norms the
density tends to zero, and the target here asks for three hits.
Scale the interval length with the mean spacing and replace the
fixed-pattern cancellation by a comparison of void probabilities
at several lengths.

Let \(b(n)\) indicate positive sums of two squares, allowing zero
coordinates, and let \(K\) be the Landau–Ramanujan constant. For fixed
\(s>0\) set \(k_X(s)=\lfloor s\sqrt{\log X}/K\rfloor\) and
\(Z_s(t)=\sum_{j=1}^{k_X(s)}b(t+j)\) for uniform integer \(t\in(X,2X]\).
Let \((z)_j=z(z-1)\cdots(z-j+1)\). The existing mean and second
factorial moment give
\[
 \mathbb E Z_s\longrightarrow s,\qquad
 \mathbb E(Z_s)_2\longrightarrow s^2.                         \tag{2}
\]
The second limit uses the stated ZF and G analytic inputs. The
fourth moment bound \(\mathbb E Z_s^4\ll_s1\) follows independently
from the upper sieve. The deduction below uses only these moment
properties and that the intervals are translates of a single set.
It neither estimates the signed centered triple sum nor assumes
any triple factorization.
[Pair moments, §7.3](CONNORS_KEATING_INDICATOR_CORRELATIONS.md#73-the-mean-second-factorial-moment-and-variance)
· [Fourth-moment bound](AVERAGED_GAUSSIAN_NORM_TRIPLES.md#a-uniform-upper-sieve-for-the-window)
· [ZF and G](PRIMITIVE_NORM_BINARY_IMPORTED_INPUTS.md)

### The all-actual fourth-moment bound

Here is the upper sieve restricted to the actual indicators used in
both lower bounds. Let \(\ell\in\{1,2,3,4\}\), put
\(H=k_X(s)\), and fix distinct shifts \(a_1,\ldots,a_\ell\in[1,H]\).
Let \(D=\prod_{i<j}|a_i-a_j|\), with empty product \(1\).
For an inert prime \(p\equiv3\bmod4\) outside the divisors of \(D\)
and with \(p>2\ell\), the \(\ell\) conditions \(v_p(t+a_i)=1\)
exclude disjoint sets of \(p-1\) residues modulo \(p^2\).
Every actual sum of two squares avoids valuation one: if \(p\mid x^2+y^2\),
the nonsquareness of \(-1\) forces \(p\mid x,y\).
The local excluded proportion is \(\ell(p-1)/p^2\).
For squarefree \(d\) the CRT gives at most \(\ell^{\omega(d)}d\)
excluded residues modulo \(d^2\) for the corresponding sieve atom.
This also bounds its interval-count remainder.

Apply the upper fundamental lemma at level \(X^{1/100}\) with
prime cutoff \(z=X^{1/(100S)}\) for a sufficiently large fixed \(S\)
depending on \(\ell\). The local proportions are at most \(\ell/p\),
so the sieve dimension is at most \(\ell\) by Mertens' bound.
The weighted remainders through that level sum to
\(O_\ell(X^{1/50}(\log X)^{O_\ell(1)})\), since the weight is
\(3^{\omega(d)}\) and a fixed factor of \(\omega(d)\) in an exponential
is bounded by a fixed-order divisor function. Comparing the retained
Euler factors with \((1-1/p)^\ell\) and restoring omitted primes gives
\[
 \sum_{X<t\le2X}\prod_{i=1}^{\ell} b(t+a_i)
 \ll_\ell
 X(\log X)^{-\ell/2}\prod_{p\mid D}(1+2\ell/p)+X^{1/40}.
\]
The product exponent follows from
\(\sum_{p\le z,\ p\equiv3\bmod4}p^{-1}
=\tfrac12\log\log z+O(1)\). For this last form,
write \(\sigma=1+1/\log z\). Replacing the partial prime sums of \(1/p\)
or \(\chi_4(p)/p\) by their full \(\sigma\)-sums costs \(O(1)\)
by the Chebyshev upper bound and partial summation at \(z\).
The Euler logarithms give \(\log\zeta(\sigma)=\log\log z+O(1)\)
and \(\log L(\sigma,\chi_4)=O(1)\), the latter also following
from the alternating series for \(L(\sigma,\chi_4)\).
[Richert, upper fundamental lemma, Theorem 11.7](sources/richert-sieve-methods-1976.pdf)

For \(\ell\ge2\) let \(v=\binom{\ell}2\) and
\(j(n)=\prod_{p\mid n}(1+2\ell/p)\). The divisor expansion of
\(j(n)^v\) has nonnegative squarefree coefficients \(c_d\),
multiplicative with \(c_p=O_\ell(p^{-1})\).
Thus \(\sum_{n\le H}j(n)^v\le H\prod_p(1+c_p/p)=O_\ell(H)\).
For each pair of shifts the sum of \(j(|a_i-a_j|)^v\) over distinct
tuples is \(O_\ell(H^\ell)\) by counting tuples with that gap.
The collision product is bounded by \(\prod_{i<j}j(|a_i-a_j|)\),
whose sum is \(O_\ell(H^\ell)\) by Hölder. The case \(\ell=1\)
has no collision product. Since \(H\ll_s\sqrt{\log X}\),
the first four factorial moments of \(Z_s\) are uniformly bounded.
The ordinary fourth moment follows by grouping repeated indices.
The first moment limit in (2) is classical: Landau's counting
asymptotic gives the mean over \((X,2X]\), and translating the
interval endpoints by at most \(H\) contributes \(O(H^2/X)\)
after summing over shifts.

**Proposition.** Under (2), for every fixed
\(\lambda>(3+\sqrt3)/6\) a positive lower proportion of the
intervals of length \(k_X(\lambda)\) contain at least three sums of
two squares. For \(\lambda>2\), the same conclusion follows from
the first moment and the upper sieve alone, independently of ZF
and G.

Write \(V_X(k)\) for the fraction of starts with no hit in the next
\(k\) integers. Periodize the sample on the circle of length \(N_X\)
formed by the consecutive integers in \((X,2X]\). For
\(k\ll\sqrt{\log X}\), periodizing changes a void probability by
\(O(k/N_X)\) and each moment of order \(j\le4\) by \(O(k^{j+1}/N_X)\).
If the circle contains hits, the void count is
\(\sum_g(g-k)_+\), where \(g\) ranges over the successive cyclic
distances between hits. Each term is convex in \(k\). With no hits
the void count is constant. Hence the periodized void probability
extends to a convex function by linear interpolation. In the
following inequalities write \(V_X(s)\) for its value at \(k_X(s)\).

Put \(Q(s)=1-s+s^2/2\) and \(M_X(\lambda)=\mathbb E(Z_\lambda)_3\).
The pointwise inclusion–exclusion inequalities
\[
 1-z+\binom z2-\binom z3\ \le\ \mathbf1_{z=0}\
       \le\ 1-z+\binom z2
\]
give \(V_X(\lambda)\ge Q(\lambda)-M_X(\lambda)/6+o(1)\) and
\(V_X(s)\le Q(s)+o(1)\). Fix
\((3+\sqrt3)/6<\lambda<1\), take \(T=3\lambda-1\) and fix
\(0<\varepsilon<\lambda\). Convexity at
\(\lambda-\varepsilon,\lambda,T\) gives
\[
 V_X(T)\ge
 1-T+\lambda T-\tfrac12\lambda^2-\tfrac12(T-\lambda)\varepsilon
 -\tfrac16\bigl(1+(T-\lambda)/\varepsilon\bigr)M_X(\lambda)+o(1).
                                                                    \tag{3}
\]
The exact secant ratio for integer lengths converges to
\((T-\lambda)/\varepsilon\). The existing fourth-moment bound keeps
\(M_X(\lambda)\) bounded, so replacing that ratio costs \(o(1)\).

For every nonnegative integer \(z\),
\(\mathbf1_{z\ge1}\ge (5z-z^2)/6\). The only positive integers where
equality holds are \(2,3\). Apply this inequality and (2) at \(T\)
to give \(V_X(T)\le 1-(4T-T^2)/6+o(1)\). Subtracting this upper
bound from the first four terms of (3) gives the margin
\[
 \frac{6\lambda^2-6\lambda+1}{6}
              -\frac{(T-\lambda)\varepsilon}{2}.
\]
Its first term is positive under the stated hypothesis and its
second term can be made smaller by choosing a fixed
\(\varepsilon>0\). It follows that \(\liminf_X M_X(\lambda)>0\).
Hölder yields
\[
 M_X(\lambda)\le
       (\mathbb E Z_\lambda^4)^{3/4}\,
       \mathbb P(Z_\lambda\ge3)^{1/4}.
\]
For an explicit form of the lower bound, put
\(L_4(\lambda)=\limsup_X\mathbb E Z_\lambda^4<\infty\) and
\(D=(6\lambda^2-6\lambda+1)/6>0\) for
\((3+\sqrt3)/6<\lambda<1\). Choose
\(\varepsilon=\min(\lambda/2,D/(T-\lambda))\).
The margin is then at least \(D/2\), so (3) and Hölder give
\[
 \liminf_X\mathbb P(Z_\lambda\ge3)\ \ge\
 \frac{1}{L_4(\lambda)^3}
 \left(\frac{3D}{1+(T-\lambda)/\varepsilon}\right)^4 .
\]
Monotonicity in the window length gives positive lower proportion
for the remaining larger \(\lambda\). When \(\lambda>2\), use
\(Z_\lambda\le2+Z_\lambda\mathbf1_{Z_\lambda\ge3}\) and
Hölder instead. The unconditional first moment and upper sieve give
\[
 \liminf_X\mathbb P(Z_\lambda\ge3)\ \ge\
 \frac{(\lambda-2)^{4/3}}{L_4(\lambda)^{1/3}} .
\]

The two-scale argument is needed for this range. For
\(0<\lambda\le1\), a variable supported on \(\{0,1,2\}\) with masses
\(1-\lambda+\lambda^2/2,\lambda-\lambda^2,\lambda^2/2\) has
precisely the limits in (2) at that single scale and has vanishing
third factorial moment. Even at \(\lambda>1\), the first two
moments do not determine the third. For example a Poisson variable
of mean \(\lambda\) can be perturbed on \(0,1,2,3\) by adding a
sufficiently small multiple of \((-1,3,-3,1)\) to its masses. This
preserves normalization and both moments in (2) but changes the
third factorial moment by six times the multiple. Likewise,
fixed-shift pair probabilities alone can miss fixed-shift triple
probabilities. Independently assign each block of three integers a
uniform even-parity binary string and choose a single uniform
translation phase for the block partition. The resulting stationary
process has density \(1/2\) and all distinct pair probabilities
\(1/4\). A consecutive triple has probability \(1/12\) rather
than the independent value \(1/8\), since one phase places it
inside a single block. These examples concern the specified
information and make no arithmetic prediction for Gaussian norms.

This establishes a lower-count subgoal of
GAP-ARITH-CK-AVERAGED-CENTERED-TRIPLES. The third factorial moment
asymptotic, its equivalent signed remainder, and a specified
fixed-gap triple asymptotic remain open. For the C12 shared-summand
problem the existing positive-proportion natural-scale theorem
already uses a proved fixed-gap pair lower bound and the independent
outer target average. The void argument requires a window process
with the moment limits (2), which the pointwise C12 target does not
supply. [C12 shared-summand theorem and precise completion target](SHARED_SUMMAND_C12_CORRELATION.md)
