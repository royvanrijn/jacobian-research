# Arithmetic completion and radical growth for least-repaired shifts

Throughout take \(n=X^L\), with \(L>0\) fixed before the
height limit, \(M=\lceil K/L\rceil\), and keep genuine
positive primitive original rows and their three radical
caps. Choose one original marked role \(T\), and a different
original role \(S\), with the fixed signs of the original
additive relation. Define
\[
 Z(T)=\prod_{p^e\parallel T}p^{(M-e)_+},\qquad
 H(T)=TZ(T),\qquad U=SZ(T).
\]
The marked role may be an original summand. Exact \(T,S\)
recover the original row. From an exact repaired letter
\(U\), first choose \(Z(T)\mid U\) at divisor cost. Its
prime powers determine the deficient primes of \(T\) and
their original valuations. The original nondeficient
support has radical at most \(X^{1/M}\le n^{1/K}\),
and after it is chosen its powers cost a height subpower.
Thus each exact \(U\) has at most \(n^{1/K}X^\varepsilon\)
original completions, uniformly at fixed \(K,L\).

## PS1. One weighted physical shift on diffuse marked inputs

Assume the marked role on the first rows has a fixed
positive amount \(\beta\) of raw mass in a compact band
of positive scaled degrees, with each prime in the
retained part of this band of arbitrarily small raw
price. A positive nonatomic component in a refined
original profile supplies this after trimming finitely
many large atom neighborhoods. Require the shared
original supports of the two rows to have raw mass
\(\le\xi\), with \(\xi<\beta/4\). For each \(\eta>0\),
sufficiently large fixed \(K\) and then sufficiently
small fixed \(L\), every prescribed exact shift
\[
                        A U_1-B U_2=F
\]
has at most \(n^\eta X^\varepsilon\) such pairs, where
\(A,B\) are prescribed nonzero integers of bounded
polynomial height, provided the retained band on row
one outside the support of \(B\) still has mass at least
\(\beta/2\). This last condition is automatic for \(B=1\).
For a prescribed coefficient pool, its failure on an
original near-full cohort costs a power gap by EB5.
A varying coefficient dictionary must still be charged.

Choose \(\mu>0\) small relative to \(\eta,\beta\). On
row one retain private prime powers from the repaired
marked band, of total raw price between \(\mu/2\) and
\(\mu\), using the small individual price condition and
omitting all shared primes and primes of \(B\).
Call their product \(D\). Fix it at \(n^\mu X^\varepsilon\)
cost. If the upper scaled edge of the original band is
\(b\), then \(D\) is \(k\)-full for
\(k=\lfloor M-b/L\rfloor\) and
\[
                 n^{k\mu/2}\le D\le n^{M\mu}.
\]
It divides \(U_1\), is prime to \(B U_2\), and therefore
also prime to \(F\); in particular zero shift gives no
solutions in this retained sector. The displayed shift
puts \(U_2\) in one fixed unit class modulo \(D\).

Write \(R_2=\prod_{p^e\parallel T_2,\ e<M}p\) and
\(W_{2,\mathrm{def}}=\prod_{p^e\parallel T_2,\ e<M}p^e\).
Then
\[
                    U_2=R_2^M (S_2/W_{2,\mathrm{def}}).
\]
The root is at most \(n^C\) for the fixed original
individual cap \(C\). The wide rational factor has
height at most \(X\). Both factors are units at \(D\).
Use UW1 with one narrow root and amplification one.
Its height constant may be \(2C/\mu\), and its wider
exponent relative to \(D\) is at most \(4/(K\mu)\) for
large \(K\). Thus the product condition of UW1 holds
for sufficiently large fixed \(K\), then \(k\) is
sufficiently large by taking small fixed \(L\).
The count costs at most \(n^{O(\mu)+1/K}X^\varepsilon\):
indeed the determinant term \(D^{B_1/k}\) is at most
\(n^{O(\mu)}\), and the exact-fibre factor
\(D^{\rho/M}\) is at most \(X^{1/M}\le n^{1/K}\)
using \(\rho=\log_D X\). Overlap of the wide factor
with the narrow root is allowed in UW1. This counts
the possible exact \(U_2\)'s. The shift then fixes
\(U_1\); recover the original rows at \(n^{2/K}\)
additional cost. Taking \(\mu\) small and \(K\) large
absorbs all these losses in \(\eta\).

The overlap restriction discards only power-smaller
pairs in an original near-full balanced cohort, by
anchoring one original row and applying EB5 to its
product support. Thus an exponent-full pair relation
with a marked nonatomic spectrum has an exact
(unweighted or fixed weighted admissible) physical
shift dictionary of exponent at least \(2/3-o(1)\).

## PS2. Bounded-radical primitive shifts are power-small

For unweighted shifts put \(g=\gcd(U_1,U_2)\) and
\(Z=|U_1-U_2|/g\). On the low-overlap pairs \(g\)
is supported on actual shared original primes, so
choosing it costs at most \(n^\xi X^\varepsilon\).
The exact primitive auxiliary \((U_1/g,U_2/g,Z)\),
in signed additive orientation, recovers the original
pair at \(n^{2/K}X^\varepsilon\) after fixing \(g\).
For a marked role with positive original nonatomic
mass, no exponent-full pair sector can have
\(\operatorname{rad}_{P^c}(Z)\le n^B\) for any fixed
\(B\), with any fixed prescribed polynomial-height
prime pool \(P\). More precisely, consider a sequence of balanced pair
strata with a common limiting original profile, marked
nonatomic mass bounded below on the first rows, and
\(K\to\infty\). At each \(K,L\) take the height limit
first, with \(L\) sufficiently small at fixed \(K\).
The restricted counts cannot approach exponent \(2/3\).
Equivalently on such profile strata there is a uniform
positive gap for sufficiently large \(K\) and sufficiently
small fixed \(L\), with thresholds and gap allowed to
depend on \(B\) and the profile margins. For the assertion
about two specified roles below one may pass to an
arbitrary limiting original profile.

Here is the height-first contradiction. Use low
overlap \(\xi\to0\), \(K\to\infty\) and \(L\to0\),
always taking the height limit at fixed choices.
Fix \(g\) and actual supported \(Z\)-powers in \(P\)
before the inner count, at the stated costs, and
take a dyadic shell \(Y\) of the auxiliary height.
The repaired operand heights have numerical
\(\log_XY/K\to1/3\): each unnormalized operand is
at least \(\operatorname{rad}(T)^M/X\) and at most
\(X\operatorname{rad}(T)^M\), and \(g\le
X n^{M\xi}\) by splitting its primes into the
original compared part and the multiplier part.
In the \(Y\)-based degree scale \(eL/\log_XY\),
the two repaired operands have purely atomic limits.
The original \(S\)-prime mass concentrates at zero
by its bounded \(X\)-height. The repaired \(T\)-prime
mass concentrates at degree \(3\), since the original
primes with valuation at least \(\delta M\) have raw
mass \(\le1/(\delta ML)\), and \((M-e)L/\log_XY\to3\)
on the remainder as \(\delta\downarrow0\). The removed
\(g\)-raw is vanishing.

The cheap full lock \(D\) on the first operand used
above has positive numerical depth also in the
\(Y\)-scale, and order \(\asymp \log_X Y/L\).
Apply the coefficient-uniform atomic-window argument
of GE1 with that lock: the second operand is wholly
atomic, and the exact supported part of \(Z\) is a
fixed coefficient. The primitive auxiliary row
population has exponent \(2/3-o(1)\), so the outside
nonatomic mass of \(Z\) is at least \(2/3\), taking
the purchase and grouping losses arbitrarily small.
But its positive nonatomic mass now supplies its own
cheap full divisor. This releases both atomic operands
by the same argument and gives a contradiction.
There is no balance hypothesis on the auxiliary
three caps in this use of GE1.

In an arbitrary balanced endpoint at least two original
roles have positive nonatomic mass. Hence any specified
two distinct marked repair roles cannot both have bounded
normalized outside radicals on an exponent-full common
pair sector. For each repaired role the compared
summand can be chosen among the other roles, holding its
position fixed. This statement follows uniformly by
taking a contrary sequence, using a common refinement of
its original pair profiles, and selecting one of the two
marked roles with positive limiting nonatomic mass.

For fixed nearby repair orders, there is also direct
fresh support control. If
\(F_m=Z_m(T_1)S_1-Z_m(T_2)S_2\) and \(m<k\),
every prime outside both original products dividing
both \(F_m,F_k\) divides
\[
                    Z_k(T_1)/Z_m(T_1)-Z_k(T_2)/Z_m(T_2).
\]
Eliminate the two \(m\)-scale unit operands to see this.
The integer on the right is nonzero on low-overlap
pairs with positive private deficient mass, and its
absolute value is at most \(n^{(k-m)C}\).
Thus a fixed number of adjacent scales can share only
a bounded raw intersection in the fresh support, while
the individual primitive shift floors above can be
forced past any fixed budget.

## PS3. Two exact physical shifts, even without a diffuse role

Write the two repaired opposite roles as \(U,V\),
with the original signed sum giving the repaired
marked value. Fix their two exact differences \(F,G\);
their corresponding signed combination is
\(\Delta=H(T_1)-H(T_2)\). On balanced pairs satisfying
the low shared-support condition, \(\Delta\ne0\)
because equal repaired marked values have the same
original target support. Put
\[
                H(T_i)=W_i R_i^M,\qquad R_i=\operatorname{rad}(T_i).
\]
The overfull coefficient product has radical at most
\(X^{2/(M+1)}\), so enumerate the exact \(W_i\) at
\(n^{2/K}X^\varepsilon\) cost. For each coefficient pair
the fixed equation
\(W_1R_1^M-W_2R_2^M=\Delta\) has subpower many
positive solutions. Divide the root gcd (its \(M\)-th
power divides \(\Delta\)) and apply Evertse's
coefficient-uniform primitive Thue–Mahler bound
for the binomial of distinct projective roots,
at the fixed \(M\ge4\). Actual original \(T_i\)
are now supported on \(R_i\) and can be recorded
at subpower cost. Finally their exact multipliers
are known. The fixed \(F\) determines the original
compared coordinates in one integral progression
with step \(Z(T_2)/\gcd(Z(T_1),Z(T_2))\) in the first
coordinate. This step is greater than \(X\) if
\(K(r_- -\xi)>2\). The full target powers on private
primes have height at most \(X\), so their repaired
multiplier has size at least \(n^{M(r_--\xi)}/X\),
where \(r_-\) is the second target's lower raw cap.
Therefore there is at most one first coordinate
in the original height range. The fixed two
physical shifts have fibre at most
\(n^{2/K}X^\varepsilon\) across the entire
original profile range.

Sources: [UW1 all signed heights](ABC_QUANTITATIVE_DEGREE_WINDOWS.md#a-product-bound-for-the-admissible-wider-height),
[GE1 atomic-window release](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees),
[Evertse, Theorem 6](https://ir.cwi.nl/pub/2463/2463D.pdf#page=14),
[mixed-power fixed-level lifting](ERDOS_943_MIXED_POWER_KERNEL_GCD_WEIGHTS.md#2-an-exact-fixed-level-lifting-problem).

## PS4. Joint factor prices on same-tail repaired target pairs

The fixed-anchor argument of RRP4 gives arbitrary radical
growth for composite repair orders. On a same-tail chart the
following joint count gives a price for each factor using both
original completions.
[Fixed-anchor radical growth RRP4](ABC_REDUCED_TARGET_RADICAL_ODDS.md#generic-normalized-differences-have-large-radical)

Fix any desired integer \(j\ge1\), and fix distinct primes
\(\ell_1,\ldots,\ell_j\ge5\). Choose the repair order
\(M\) to be a multiple of their product with \(ML=K+o(1)\)
at small fixed \(L\), for example by rounding \(K/L\)
up to the next multiple. This does not change the
preceding support or height bounds at leading order.
Fix one exact overfull tail \(W_0\) of the marked role
on a near-full original cohort at cost at most
\(1/K+o(1)\) in the exponent. This is possible by
the overfull radical bound and supported enumeration.
For pairs from this cohort with small original shared
radical, write
\[
    R_i=\operatorname{rad}(T_i),\quad
    g=\gcd(R_1,R_2),\quad r_i=R_i/g,\qquad
    \frac{H(T_1)-H(T_2)}{\gcd(H(T_1),H(T_2))}
                                =r_1^M-r_2^M .
\]
The roots \(r_1,r_2\) are coprime and distinct on
low-overlap pairs. At such a pair, for each \(j'\)
the homogeneous cyclotomic value
\(\Phi_{\ell_{j'}}(r_1,r_2)\) is a nonzero factor of
the primitive target difference. Its form has
degree \(\ell_{j'}-1\ge4\) and distinct complex
projective roots. A fixed exact value recovers its
primitive root pair at subpower cost by Evertse,
uniformly in the value, at fixed \(\ell_{j'}\).
Its possible values with radical outside a
prescribed polynomial-height pool at most \(n^\beta\)
therefore determine at most \(n^\beta X^\varepsilon\)
primitive root pairs by supported enumeration.
Recovering \(g\) from shared original supports costs
at most \(n^\xi X^\varepsilon\). Exact original
\(T_i\) are then supported on their known \(R_i\)
and cost only subpower enumeration. Their two
original row completions cost at most \(n^{2u+o(1)}\)
on a regular original cohort with at most
\(n^{u+o(1)}\) rows over each marked complete letter. All three original caps are retained in this
last count.

Consequently each small-radical factor sector has
at most \(n^{\beta+\xi+2u+o(1)}\) original pairs
inside the chosen same-tail chart. The chart has
pair exponent at least \(2/3-2/K-o(1)\), with
the original low-overlap exception discarded by
EB5. In particular for any fixed small
\(\delta>0\), choose \(\xi\) small and \(K\) large:
each factor has outside radical greater than
\(n^{2\alpha-\delta}\), \(\alpha=1/3-u\), on all
but a power-smaller fraction of pairs in the
chart. The omission of low-overlap exceptions
can also be made before taking the large chart:
choose \(K\) large compared with \(1/\xi\),
since the global EB5 anchored bound has a gap
of order \(\xi\). Only finitely many factors
are used.

Enlarge the prescribed pool by the nonzero pairwise
resultants of the chosen homogeneous cyclotomic
forms. Their common prime supports on the reduced
primitive roots are now disjoint: a common prime
would give a common projective root of two forms.
Hence on the retained pair population the fresh
radical of the primitive repaired target shift is
at least \(n^{j(2\alpha-\delta)}\). Since
\(\alpha\ge1/6\) at every balanced endpoint and
\(j\) is arbitrary and fixed before the height
limit, this gives arbitrarily large fixed raw
floors without any nonatomic hypothesis. One
may also omit both original product supports
at the expense of their two original caps, whose
combined raw price is bounded. The conclusion concerns the retained same-tail
population at the specified composite orders.

Sources: [Evertse, Theorem 6](https://ir.cwi.nl/pub/2463/2463D.pdf#page=14)
and [primitive cyclotomic factor recovery](ABC_REPAIRED_WEIGHTED_TARGET_DIFFERENCES.md#several-nonlinear-factors-at-one-cost).

