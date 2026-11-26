# Polynomial relations among least-repaired target heights

Fix \(L>0\), let \(n=X^L\), and take the height limit at
this fixed budget. For fixed \(K>0\) put \(M=\lceil K/L\rceil\).
The least target repair of any \(c\le X\) has height
\[
 H(c)=\prod_{p\mid c}p^{\max\{M,v_p(c)\}}
      =W(c)\operatorname{rad}(c)^M,
 \qquad W(c)=\prod_{p:\,v_p(c)>M}p^{v_p(c)-M}.
\]
Only the tail \(W(c)\) needs a coefficient record.
A bounded-order polynomial relation then counts the
prime packets it can see at nearly subpower cost in \(n\).
The remaining packets act by coordinatewise dilation.
This statement concerns the actual least repair of each
target and the corresponding primitive quotient.
[Least one-target repair (MPR)](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)

## Prime packets invisible to a polynomial relation

Fix \(k\ge2\), \(B>0\), and a finite support of distinct
monomials in \(k\) variables with at least two monomials.
Let \(G(\mathbf t)=\sum_{\nu=1}^s A_\nu\mathbf t^{\mathbf e_\nu}\)
have that support, with nonzero rational coefficients.
An equality \(G(\mathbf t)=0\) is nondegenerate if no
nonempty proper subsum vanishes on \(\mathbf t\).
Call a nonempty set of indices \(S\subseteq\{1,\ldots,k\}\)
invisible if \(\sum_{i\in S}e_{\nu i}\) is independent of
\(\nu\). Define \(V_G\subseteq\mathbb Q^k\) to be the span
of the indicator vectors of the invisible sets and let
\(r=\dim V_G\). This depends only on the monomial support.

For every \(\eta>0\) and all sufficiently small fixed
\(L>0\), the nondegenerate solutions with
\(t_i=H(c_i)\), \(c_i\le X\), and
\(\operatorname{rad}(c_i)\le n^B\), lie in at most
\[
                          O(n^{k/K+\eta}X^\varepsilon)       \tag{HP1}
\]
classes under the following equivalence. Two positive
integer vectors are equivalent if, at every rational
prime \(p\), the difference of their \(p\)-valuation
vectors belongs to \(V_G\). The bound is uniform in
the nonzero coefficients. The threshold on \(L\)
depends on \(k,B,K,\eta\) and the monomial support.
As usual the \(\varepsilon>0\) allowance can be taken
arbitrarily small at fixed \(L\). In particular, if all
\(t_i\) belong to a dictionary \(\mathcal H\) of size \(D\),
the tuple count is
\[
                            O(D^r n^{k/K+\eta}X^\varepsilon).
                                                               \tag{HP2}
\]
If the slots have different dictionaries, the factor \(D^r\)
can be replaced by the product of the sizes in any \(r\)
coordinate slots whose projection is injective on \(V_G\).

To prove the result, write \(W_i=W(c_i)\) and
\(R_i=\operatorname{rad}(c_i)\).
Since \(W_i\le X\) and
\(\operatorname{rad}(W_i)^{M+1}\le c_i\), the \(k\)
tails cost at most \(X^{k/(M+1)+\varepsilon}\) choices.
This uses radical enumeration and the fixed-support Rankin
estimate, uniform for prescribed prime products of polynomial
height. The exponent \(k/(M+1)\) is at most \(kL/K\).

For each nonempty \(S\subseteq\{1,\ldots,k\}\), let \(x_S\)
be the product of primes dividing exactly those \(R_i\)
with \(i\in S\). These variables are squarefree, pairwise
coprime, and at most \(n^B\). Substitute
\[
 t_i=W_i\prod_{S\ni i}x_S^M
\]
in the relation. Every invisible variable factors out
of the full relation and all its subsums with the same
exponent. Each remaining variable has exponent spread
at least \(M\) among the monomials. For \(s\ge3\), equation (EB1) for
nondegenerate power relations counts the remaining
pairwise coprime variables. There are at most \(2^k-1\)
of them, their exponents are fixed at the chosen budget,
and the bound is uniform in the tail coefficients.
For at least two remaining variables the count at fixed
tails is at most
\[
          O\!\left((n^B)^{\,C_{k,s}M^{-1/(2(2^k-2))}
                                                   +\varepsilon'}\right).
\]
If no variable remains there is a single empty record.
With one variable, apply the univariate version. Two
different original exponent vectors cannot yield the
same exponents at every visible variable, since the
singleton of an index where they differ is visible.
Thus the modified relation is nonzero as a polynomial.
In all cases the count at fixed tails is
\(O(n^\eta X^{\varepsilon})\) for the stated small
fixed budgets, with the allowances split between the
steps. This proves (HP1): after the tails and visible
variables are fixed, variation of the others changes
each prime valuation vector only in \(V_G\).
[Coefficient-uniform nondegenerate power release (EB1)](ABC_NONDEGENERATE_POWER_COLLISIONS.md)

For two monomials a direct count applies. After factoring
out the invisible variables, write the equation as
\(\prod_S x_S^{M d_S}=C\) over the visible sets, where the
fixed integers \(d_S\) are nonzero. A solution forces the
fixed rational \(C\) to have polynomial height by the
root caps, regardless of the original coefficient heights.
The roots \(x_S\) are squarefree and pairwise coprime,
so their union is exactly the prime support of \(C\).
Each such prime can be assigned to at most \(2^k-1\)
packets with a compatible exponent. This gives
\(X^\varepsilon\) records by the uniform bound on the
number of distinct prime factors at polynomial height.

Choose \(r\) coordinate positions whose projection is
injective on \(V_G\), by taking independent rows of a
basis matrix. Inside any one equivalence class, fixing
the actual heights in those positions fixes the whole
vector, since their valuation differences are then
zero at all primes. This proves (HP2). The argument
applies to \(r=0\) with an empty coordinate record.

## Shifted sums and true additive diagonals

For \(G=\sum_i A_i t_i\) with all \(A_i\ne0\),
the only invisible nonempty set is the full index set.
A homogeneous nondegenerate linear relation thus has
at most \(O(n^{k/K+\eta}X^\varepsilon)\) rational
projective directions and at most
\[
                   O(D n^{k/K+\eta}X^\varepsilon)             \tag{HP3}
\]
tuples from a dictionary of size \(D\). Indeed, \(V_G\)
is the span of the all-ones vector, and two integer
height vectors in the same class differ by a common
positive rational scalar. For a fixed nonzero shift
\(N\), take \(G=\sum_i A_i t_i-N\). All nonempty
packets are visible and (HP1) now counts actual
height tuples. Both conclusions are uniform in the
coefficients and, in the shifted case, in \(N\).
For one slot a fixed nonzero shifted relation determines
its height directly.

Suppose instead the nonzero shift can vary, the rational
weights have polynomial height, and the reduced numerator
of \(N=\sum_i A_i H(c_i)\) has radical at most \(n^b\)
outside a prescribed polynomial-height prime pool.
Then there are at most \(O(n^{b+k/K+\eta}X^\varepsilon)\)
height tuples for which the shifted relation is
nondegenerate. The denominator primes of \(N\) belong
to the fixed product of the weight denominators. The
numerator and denominator have polynomial height,
since \(R_i^M\le n^{BM}\) and \(LM\le K+L\).
The number of rational \(N\)'s is thus at most
\(n^b X^\varepsilon\), by the fixed-support count.
In particular the nondegeneracy assumption is automatic
for positive weights.

For each fixed \(h\ge2\) consider the ordered \(2h\)-fold
energy of \(\mathcal H\), the number of equal sums of
\(h\) dictionary values on each side. Every relation
partitions into minimal nonempty vanishing blocks. A
block of size two matches a value on one side with
the same value on the other. Therefore a tuple outside
the permutation diagonal has a block of size at least
three and at most \(h-1\) blocks in total. Apply (HP3)
to each block of size at least three, and choose one
height per block of size two. The number of partitions
is bounded for fixed \(h\). This gives the bound
\[
 \#\{\text{off-diagonal equal sums}\}
             \ll D^{h-1} n^{2h/K+\eta}X^\varepsilon .        \tag{HP4}
\]
For a fixed \(u>0\), if \(D\ge n^u\), one can choose
\(K\) large in terms of \(h,u\) and then \(\eta\) small.
Taking the height limit at the stated fixed budget
gives a power loss from \(D^h\). The permutation
diagonal has size \(h!D^h+O_h(D^{h-1})\), because
all \(h\) values are distinct except on \(O_h(D^{h-1})\)
ordered lists. Consequently the \(h\)-fold sumset has
cardinality \((1/h!+o(1))D^h\) by Cauchy–Schwarz and
the count of unordered lists. Every near-minimum
target dictionary at the balanced endpoint is covered
by choosing \(K\) large enough. This uses no
regularity of its exact valuation distribution.

In particular, when \(D\ge n^u\) with fixed \(u>0\),
a fixed nonzero difference of two equal-length \(h\)-sums
has at most \(O(D^{h-1} n^{2h/K+\eta}X^\varepsilon)\)
representations, including degenerate ones.
Partition the \(2h\) signed heights and the fixed
nonzero difference into minimal vanishing blocks. The
block containing the constant uses at least one height
and costs no factor of \(D\). At most \(h-1\) homogeneous
blocks remain, and the shifted and homogeneous bounds
multiply with total tail cost at most \(2h/K\).

There is also a bound for weighted target-difference
dictionaries. Take positive rational weights \(A,B\) of
polynomial height and any prescribed polynomial-height
prime pool \(P\). For any size-\(D\) dictionary of repaired
heights the number of pairs with \(A H_1-B H_2=0\), or
with a nonzero reduced numerator of outside-\(P\)
radical at most \(n^b\), is
\[
                   O(D+n^{b+2/K+\eta}X^\varepsilon).         \tag{HP6}
\]
The zero difference gives at most one second height
for each first height. The nonzero difference gives
a nondegenerate three-term equation, so the moving
shift bound applies. If \(D\ge n^d\) with fixed \(d>0\),
the nonzero outside radical generically exceeds
\(n^{2d-2/K-\rho}\) for every fixed \(\rho>0\), in the
sense of a power-small exceptional set of pairs. Take
\(\eta\) and \(\varepsilon/L\) smaller than \(\rho\).
An allowed dictionary of \(n^\gamma\) pairs of weights
reduces the threshold by \(\gamma\) for the nonzero
differences, by summing over the dictionary. The
same summation contributes \(O(n^\gamma D)\) on
the zero locus. Thus a large repaired dictionary by
itself supplies few cheap weighted differences of
the type used in kernel recovery.
[Repaired weighted target differences](ABC_REPAIRED_WEIGHTED_TARGET_DIFFERENCES.md)

## Transfer back to primitive abc rows

Assume now \(\mathcal E\) is a family of positive primitive
\(a+b=c\) with \(X/2<c\le X\), and each original
complete letter has radical at most \(n^{1/3+\kappa}\).
Write \(T=|\mathcal E|\) and let \(F\) be the largest
number of primitive rows with a given \(H(c)\).
At fixed \(H(c)\), \(R=\operatorname{rad}(H(c))\) is known
and \(c\) has exactly this support. The fixed-support
count followed by the complete-letter bound gives
\[
                 F\ll n^{1/6+\kappa/2+\eta}X^\varepsilon .
\]
This argument concerns representations with the unique
least multiplier. It introduces no arbitrary common
multiples. The repaired-height dictionary therefore
has size at least \(T/F\).
[Fixed-letter count and repaired target radical fibres](ABC_FIXED_LETTER_RADICAL_CAPS.md#size-of-a-minimally-repaired-target-radical-fibre)

In (HP2) choose the \(r\) projected slots as primitive
rows instead of as heights. Each surviving equivalence
class then fixes the other \(k-r\) heights, which lift
at most \(F^{k-r}\) ways. Therefore the number of original
row tuples whose repaired heights satisfy a given
nondegenerate relation is
\[
                  O(T^r F^{k-r} n^{k/K+\eta}X^\varepsilon).
                                                               \tag{HP5}
\]
For \(r<k\) this gives a uniform positive loss from
the full tuple exponent \(k/3\) near saturation, after
choosing \(K\) sufficiently large, then the cap and
count tolerances sufficiently small. It works throughout
the balanced endpoint. The same conclusion applies
if the prescribed coefficient dictionary has size at
most \(n^\gamma\), with a loss reduced by at most
\(\gamma\). For nondegenerate homogeneous linear
relations the loss can approach \((k-1)/6-k/K\).
From (HP4), a family of saturated size also has almost
no genuinely off-diagonal fixed-order equal sums on
its repaired target dictionary, in the sense of the
power-saving error in (HP4). Equal repaired heights can
still have several original representations.


There is a coefficient entropy consequence of (HP5).
For any deterministic assignment of a coefficient
record with that monomial support to tuple data,
suppose it gives a nondegenerate relation on a set
\(\mathcal A\subseteq\mathcal E^k\). On the uniform law
of \(\mathcal A\), the entropy of the coefficient record,
even when records are taken only up to nonzero rational
scaling, is at least
\[
       \log |\mathcal A|-r\log T-(k-r)\log F
                -(k/K+\eta)\log n-O(\varepsilon\log X+1).
\]
Indeed the largest coefficient fibre is bounded by
(HP5) for a fixed representative of the projective
coefficient vector. The bound also applies to the
negative logarithm of the largest record probability.
Near saturation, any such rule valid on a full-exponent
tuple population uses at least
\((k-r)/6-k/K-o(1)\) normalized coefficient entropy
after the cap and bound allowances are made small.
A nondegenerate relation with a cheap coefficient
dictionary would therefore need special arithmetic
origin in the original additive graph.

The bounds address an exact arithmetic condition on
tuples of repaired heights. To turn (HP5) into a
one-row count, an extraction from an actual near-full
\(a+b=c\) family would have to produce a sufficiently
large population of repaired-height relations with
\(r<k\) and affordable coefficients, or impose an
affordable nonzero shifted output cap. Additive
relations on the original complete labels need a
separate argument to survive the row-dependent repair.
[Minimal-repair count equivalence](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)
