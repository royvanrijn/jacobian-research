# Large digit contrast images under a small additive congruence

In the balanced setting, a small modulus can support recovery
of almost all the nonatomic radical of two letters. To bound
the number of these records, the original sum at full precision
is important. A family of actual integers can satisfy all
three radical caps, have a prescribed full divisor on the
third letter, satisfy the sum modulo its small divisor, and
come within a relatively short real interval of the sum. The
small-modulus contrast image can still be much larger than
the one-third target.

## The local model

Fix \(0<\delta<1/100\), let \(w=3\delta\), and set \(n=X^L\).
The height limit is always taken at fixed \(L>0\) before
\(L\) goes to zero. Fix a small \(0<\tau<1/4\).

**SMO1.** There are families of positive, pairwise coprime
integer rows \((U,V,T)\), with each coordinate comparable to
\(X\) and at most \(X\), such that
\[
 \operatorname{rad}(U),\operatorname{rad}(V),\operatorname{rad}(T)
 \le n^{1/3}.
\]
The normalized exact-degree measures of \(V,T\) converge to
nonatomic measures of mass \(1/3\). The \(U\)-measure converges
on a subsequence of budgets to a nonatomic part of mass
\(1/3-w\) and one positive atom of mass \(w\). There is a common
full divisor \(D'\) on \(V\) from compact positive valuation
levels, with normalized radical weight \(\tau/3+o(1)\) and
positive limiting numerical weight. The family has a common
divisor \(D\mid D'\) with \(\log_n D=\delta+o(1)\). Its rows obey
\[
           U+V-T\equiv0\pmod D,\qquad
           |U+V-T|\ll_L Xn^{-\delta/4}.
                                                               \tag{SMO1}
\]
The small-budget and inner height errors in the weight claims
follow the height-first convention. One can hold \(V\) fixed
or allow it a conditional dictionary of at least
\(n^{(1-\tau)/3-o_X(1)-o(1)}\) values per retained choice of the
two diffuse cores. The rows can be chosen to fail the corresponding congruence
modulo \(D'\).

For the construction take \(m=\lfloor L^{-1/2}\rfloor\),
an odd prime \(d\) between \(2/L\) and \(4/L+2\), and
distinct integers \(e_1,\ldots,e_m\) with \(Le_j\) at
distance \(O(L)\) from \(2+2(j-1/2)/m\). Bertrand's
postulate gives \(d\). An adjustment of a single \(e_j\)
by one ensures \(d\notin\{e_j\}\) since the consecutive
gaps tend to infinity. Put
\[
 t_j=Le_j,\quad \overline t=m^{-1}\sum_j t_j,\quad
 s^2=m^{-1}\sum_j(t_j-\overline t)^2 .
\]
Then \(\overline t\to3\) and \(s^2\to1/3\). For a desired
mass \(\alpha\) and numerical weight \(\beta\), give slot
\(j\) normalized radical weight
\[
    a_j(\alpha,\beta)=\frac{\alpha}{m}
        \left(1+\frac{\beta/\alpha-\overline t}{s^2}
                         (t_j-\overline t)\right).
\]
Its masses sum to \(\alpha\) and its numerical weights
sum to \(\beta\). For \(V,T\) use \((\alpha,\beta)=(1/3,1)\);
for the core \(U_0\) of \(U\) use
\((\alpha,\beta)=(1/3-w,1-Ldw)\).
The weights are positive and \(O(1/m)\) for every sufficiently
small \(L\). Indeed the desired core mean differs from three by
\(w|3-Ld|/(1/3-w)\), at most \(0.11\) at late budgets.
Each slot consists of a prime of size a fixed multiple of
\(n^{a_j(\alpha,\beta)}\) raised to \(e_j\), with the prime
chosen from a fixed positive-fraction subinterval. The
prime number theorem supplies \(n^{a_j-o_X(1)}\) choices in
each such interval. All earlier selected primes can be
excluded at negligible cost. The ordinary prime number
theorem and its fixed-progression version are used only
after fixing \(L,d,m\).

For \(T\) choose the constants of these prime intervals so that
\(3X/4\le T\le7X/8\), for \(V\) so that
\(X/4\le V\le X/3\), and for \(U_0\) so that
\[
        (3/4) X^{1-Ldw}\le U_0\le X^{1-Ldw}.
\]
For example if \(E=\sum_j e_j\), for an interval
\([aX^\beta,bX^\beta]\) use the bounds \(a^{1/E}\) and
\(b^{1/E}\) on the ratio of each base to its nominal
\(n^{a_j}\). All prime bases are at most their nominal
powers of \(n\).

Fix the primes of \(V\) in the first \(\lfloor\tau m\rfloor\)
slots in residue \(2\bmod d\), using the prime number theorem
in this fixed progression. Their full powers form \(D'\). Its
normalized radical weight tends to \(\tau/3\), and its numerical
weight tends to \((2\tau+\tau^2)/3\). Select \(D\mid D'\) by
multiplying individual prime factors of \(D'\), up to their
available exponents, until first exceeding \(n^\delta\).
This is possible for small \(L\). The overshoot costs
\(O(1/m)+o_X(1)\) in its base-\(n\) logarithm. Since the
primes of \(D\) are different from \(d\) and are \(2\bmod d\),
\(d\) is coprime to \(\varphi(D)\). Thus raising to the \(d\)-th
power permutes the unit classes modulo \(D\).

Select the two cores \(U_0,T\) after excluding all primes of
\(D'\). There are at least \(n^{2/3-w-o_X(1)}\) distinct
ordered core pairs, by the prime counts and the distinct
orders on each rail. Next the remaining slots of \(V\) give
a conditional dictionary of size
\(n^{(1-\tau)/3-o(1)-o_X(1)}\) after forbidding primes of
\(U_0,T\). One may instead fix all of \(V\) before selecting
the cores and forbid all its primes there. All three chosen
integers are pairwise coprime. Write \(G=U_0VT\).

Let \(b=( (T-V)/U_0 )^{1/d}\). The chosen heights ensure
\(b\asymp_L n^w\) and \(b\le(5/6)^{1/d}n^w\). There is a
unique unit class \(a\bmod D\) with
\[
                           a^d\equiv T U_0^{-1}\pmod D.
\]
Choose a squarefree integer \(W\equiv a\bmod D\) in
\([b,b+H]\), where \(H=n^{11\delta/4}\), coprime to \(G\).
Existence uniformly in the core choices follows from the next
estimate. Set \(U=U_0W^d\). Then \(U\equiv T\bmod D\), while
\(D\mid V\). For large \(X\) at fixed \(L\), \(W\le n^w\)
and \(U\le X\). The mean-value theorem gives
\[
                  |U+V-T|\ll_L XH/n^w=Xn^{-\delta/4}.
\]
Because \(W\) is squarefree and coprime to the earlier factors,
it contributes just the prescribed atom and the total radical
cap of every rail follows from its weights. The diffuse
limiting measures follow from the mesh of the \(t_j\), the
uniform \(O(1/m)\) bound on the weights, and the limiting tilt
formula along a subsequence on which \(Ld\) converges.
[Prime number theorem in fixed progressions, de la Vallée Poussin](https://books.google.com/books/about/Recherches_analytiques_sur_la_th%C3%A9orie_d.html?id=7e0GAAAAYAAJ)

### Squarefree drivers with complete coprimality

For \(B\asymp_L n^w\), \(H\le B\), a unit class \(a\bmod D\),
and \(G\) of polynomial height, the number of squarefree
\(x\in[B,B+H]\) with \(x\equiv a\bmod D\), \((x,G)=1\)
has the expression
\[
 \frac{H}{D}
  \prod_{\substack{p\mid G\\p\nmid D}}\left(1-\frac1p\right)
  \prod_{p\nmid DG}\left(1-\frac1{p^2}\right)
  +O\left(\tau_0(G)\sqrt B+\frac{H}{D\sqrt B}\right),
\]
where \(\tau_0\) is the divisor function. Indeed use
inclusion-exclusion over divisors of the part of \(G\)
outside \(D\), and insert \(\sum_{k^2\mid x,\ (k,DG)=1}
\mu(k)\) for squarefreeness. The inner sum may be
truncated at \(\sqrt{B+H}\). Each CRT count has bounded
error, and the Euler product tail costs \(O(B^{-1/2})\).
The first product is at least \(2^{-\omega(G)}=n^{-o_X(1)}\)
uniformly in the polynomial-height range, and the second
is bounded below by \(1/\zeta(2)\). In the construction
\(\log_n D\le\delta+O(1/m)+o_X(1)\), so \(H/D\) has
exponent \(7\delta/4-O(1/m)-o_X(1)\). This dominates the
squarefree sieve error of exponent \(3\delta/2+o_X(1)\)
at late budgets. In particular the count tends to infinity
at each such fixed budget. The \(d\)-th power also permutes
the units modulo \(D'\), and \(D'\) is larger than \(H\) at
small \(L\). Thus at most one driver in the interval satisfies
\(U\equiv T\bmod D'\), and one can delete it. This also
removes any exact solution of \(U+V=T\).

## The contrast image

Use modulus \(D\) and no initial outside pool. A mask here
can be any prescribed exact rational product of the form
\[
       \prod_{p\mid U}p^{f_U(v_p(U))}
           \prod_{p\mid T}p^{-f_T(v_p(T))},
\]
where each \(f_U,f_T\) is integer-valued and fixed at the
chosen budget. This includes the exact rational contrasts
of (MR4) and integral products of projected ratios.
At late budgets
take a fixed loss \(\zeta>0\) and any list of those contrasts
such that each selected mask has numerator times denominator
at most \(n^{\delta/3}\) and their union detects normalized
diffuse radical weight at least \(2/3-w-\zeta\) on the two
rails. The mixed-support probe theorem provides such a bounded
list with strict margins, by using a smaller precision parameter
in its admissibility condition. All words about the limiting
profile and finite root slots in that theorem apply to these
literal prime powers and squarefree drivers. Its proof concerns
the valuation measures and coprimality of the two rails and
does not impose exact additivity. One can also refine each
of the base intervals into dyadic boxes at fixed-budget cost.

**SMO2.** The image of these contrasts modulo \(D\) has size
at least
\[
                           n^{2/3-w-\zeta-o_X(1)}.             \tag{SMO2}
\]
It has the same lower bound with \(V\) held fixed. In particular
the exponent exceeds \(1/3\) whenever \(w+\zeta<1/3\).

An active driver would contribute its full radical
\(n^{w+o_X(1)}\) to the product of numerator and denominator.
The driver has disjoint prime support from both cores, so
such a contribution cannot cancel in a contrast. The
admissible masks are therefore independent of the driver
and of \(V\). Each exact rational mask is known from
its residue, since the determinant of two candidates has
absolute value at most \(n^{2\delta/3}<D\). Consequently
a fixed transcript supplies a fixed prime set detecting
at least \(2/3-w-\zeta\) of the core radical. Record the
remaining core support at cost at most \(n^{\zeta+o_X(1)}\),
using its squarefree product, and assign its primes and
the known primes to core slots at subpower cost. This
counts the ordered cores in that transcript. Dividing
the number of core pairs by this bound gives (SMO2).
[Bounded mixed-support probe theorem (MR4)](ABC_MULTI_SCALE_MIXED_RESIDUE_RECOVERY.md#a-bounded-number-of-small-modulus-probes)

Thus a uniform compression of the joint contrast dictionary
cannot follow from the three radical caps, the balanced
positive degree profiles, the existence of a common
cheap full divisor, the additive relation reduced to the
small modulus, and the displayed amount of real precision.
A count also using \(U+V\equiv T\bmod D'\), greater real
precision, or the exact integer sum has a stronger hypothesis.

## Simultaneous small congruences on all three rails

The same obstruction can impose the sum modulo one small
divisor on each original letter. Fix \(0<\delta<1/200\),
put \(w=6\delta\), and use the preceding exponent mesh and
driver order \(d\). Each rail has a core of mass
\(\alpha=1/3-w\) and numerical weight \(1-Ldw\), obtained
from the same weight formula. Each core can be chosen in
\([(9/10)X^{1-Ldw},X^{1-Ldw}]\) by fixing the constant-factor
prime intervals in advance. Fix disjoint primes in
residue \(2\bmod d\) in the first \(\lfloor\tau m\rfloor\)
slots of each core. Their full powers give common locks
\(D'_U,D'_V,D'_T\). Their normalized radical weight
\(\mu_L\) per rail is \(O(\tau)\) and bounded below by a
positive constant depending on \(\tau\). Their numerical
weights are bounded below as well. In each lock select the
small divisor \(D_U,D_V,D_T\) by first crossing \(n^\delta\)
as above. They are pairwise coprime and each has base-\(n\)
logarithm \(\delta+O(1/m)+o_X(1)\). Choose the remaining
core primes freely, with all primes distinct.

**SMO3.** The three letters admit squarefree drivers, each
of mass \(w+o_X(1)\) at order \(d\), so that the resulting
positive pairwise coprime integers satisfy
\[
   U+V\equiv T\pmod{D_U D_V D_T},\qquad
   |U+V-T|\ll_L Xn^{-3\delta/4}.
\]
The radical caps remain \(n^{1/3}\) per letter, and the
letters are comparable to \(X\) and at most \(X\).
On any selected pair of rails, using the modulus on the
opposite rail, a cheap contrast list detecting all but
\(\zeta\) of that pair's diffuse mass has image of size at
least \(n^{2(1/3-w-\mu_L)-\zeta-o_X(1)}\). This exceeds
\(n^{1/3}\) when \(\delta,\tau,\zeta\) are sufficiently small.
The bounded probe theorem again supplies such lists.

Write \(V=V_0W_V^d\) and choose its driver to give
\(V\equiv1\bmod D_U\), \(V\equiv-1\bmod D_T\). Choose
\(T=T_0W_T^d\) with \(T\equiv1\bmod D_U D_V\). These
are prescriptions of one unit class for each driver,
using the \(d\)-th power permutation and the CRT.
Use for \(W_V\) the interval
\([(1/4)^{1/d}n^w,(1/3)^{1/d}n^w]\), and for \(W_T\)
the interval \([(3/4)^{1/d}n^w,(7/8)^{1/d}n^w]\).
The squarefree progression estimate applies at a modulus
of exponent \(2\delta+O(1/m)+o_X(1)\), since \(w>4\delta\).
Each driver can be taken coprime to all cores and all
previous drivers. Thus \(T-V\) lies in a fixed positive
interval of multiples of \(X\), and is at most \(0.65X\).

Choose \(U=U_0W_U^d\) with \(U\equiv1\bmod D_V D_T\),
and \(W_U\) within \(H=n^{21\delta/4}\) above
\(((T-V)/U_0)^{1/d}\). The squarefree progression estimate
again applies because the main term exponent
\(13\delta/4-O(1/m)-o_X(1)\) exceeds the error exponent
\(3\delta+o_X(1)\). The center stays below a fixed constant
less than one times \(n^w\). This proves the real bound
and the height and radical caps. The three congruences
follow from \(D_U\mid U_0\), \(D_V\mid V_0\), and
\(D_T\mid T_0\) and the chosen unit classes. For any fixed
cores and earlier choices, one can also delete the at
most one \(U\)-driver producing exact equality.

For a chosen complementary pair, the independently
assigned varying core primes give \(n^{2(\alpha-\mu_L)-
o_X(1)}\) core pairs, even with the opposite core fixed.
Admit masks with numerator times denominator at most
\(n^{\delta/3}\). They cannot hit either driver, since
each has radical \(n^{w+o_X(1)}\) on a disjoint support.
Recover their exact rationals from the residues and
enumerate the remaining core support as in (SMO2).
This gives the stated image lower bound. Requiring
all three small congruences and using all three lists
together cannot decrease the lower bound from one list.
