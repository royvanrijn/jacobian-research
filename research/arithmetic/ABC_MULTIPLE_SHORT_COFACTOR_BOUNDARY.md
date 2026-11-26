# Boundary release with several short cofactors

A bounded set of high-order roots can be released when
each letter also has its own small integral cofactor.
The boundary argument controls the high-root supports.
After those roots are fixed, a large coefficient in the
remaining linear equation determines the small cofactors
by rational reconstruction. This replaces both real-gap
terminals in the singleton dictionary proof.

## Releasing the high roots and all three cofactors

**MCB1.** Fix \(m\ge0\), \(B>0\), and let \(n=X^L\)
at a fixed budget \(L>0\). Partition
\(m\) original roots \(x_i\) among three letters and
write
\[
 T_j=A_j z_j\prod_{i\in I_j}x_i^{e_i},\qquad
 j=1,2,3,
\]
where the positive \(A_j\) are fixed, each \(x_i\) and
\(z_j\) is a positive integer at most \(H_*=n^B\).
Require coprimality of variable coordinates across
different letters. Coordinates on the same letter
may share primes. Assume each \(e_i\ge E\), and
fix the orders before the height limit. Count actual
primitive solutions of \(T_1+T_2=T_3\) with
\(\max_j T_j>2H_*^3\). Their number is
\[
 O_{m,B,\mathbf e,L,\varepsilon}
    \left(n^{B C_m E^{-1/(2(m+2))}}X^\varepsilon\right)
 \qquad(m\ge1,\ \varepsilon>0).
 \tag{MCB1}
\]
The constant is uniform in the \(A_j\). For \(m=0\)
there is at most one solution. An absent or fixed
cofactor can be represented by \(z_j=1\) after
absorbing its value into \(A_j\).

Use the geometric dimension descent from (EB1) on
the \(m+3\) variable coordinates with height bound
\(H_*=n^B\). The equation is a nonzero polynomial,
since the factors \(z_j\) are separate variables.
The sample is nondegenerate because all three
complete terms are positive. Consider a geometric
integral component of degree \(D\) during the
descent. If some high coordinate \(x_i\) has an
affine zero on the component at which every coordinate
on the other two letters is nonzero, the
Brownawell–Masser argument of (EB1) gives
\(D\gg_m\sqrt E\). Indeed a general linear curve
section meets that boundary divisor. The letter
containing \(x_i\) vanishes to order at least \(e_i\).
Its same-letter coordinates are regular at this
affine point and may also vanish. The two other
terms are units there. The valuation of the ratio
of the first term to either of the other two is
therefore at least \(e_i\), while its function
degree on the curve is \(O_m(D^2)\). The cofactors
are coordinate functions on the curve and contribute
only \(O_m(D)\) zeros and poles.

Use the high-degree covering and proper dimension
descent of (EB1) whenever \(D\ge c_m\sqrt E\).
In dimension \(r\le m+2\) that cover has cost
\(H_*^{(r+1)D^{-1/r}+\varepsilon'}\), with the
degree constants fixed before the height limit.
Nonrational geometric components descend on
conjugate intersections as in (EB1).
On a branch with \(D<c_m\sqrt E\), every high
coordinate lacks the indicated affine boundary.
The same bounded-height interpolation and arithmetic
Nullstellensatz proof used for (SHF5) supplies,
for each \(x_i\), a nonzero integer of height
\(H_*^{O_{m,D}(1)}\) supporting all its sample
prime factors. Here the system for that coordinate
uses \(x_i=0\) and the inverse of the product of
all variable coordinates on the other two letters.
Cross-letter coprimality gives the required
inverse at every prime divisor of \(x_i\).
If the interpolation replacement lowers dimension,
descend first. Otherwise this fixes all the high
roots at subpower cost by uniform fixed-support
enumeration. This argument places no support
restriction on the cofactors.
[Exclusive-boundary argument and high-degree covers (EB1)](ABC_NONDEGENERATE_POWER_COLLISIONS.md#exclusive-boundary-release-for-nondegenerate-power-collisions)
· [One-cofactor version (SHF5)](ABC_FIXED_SHIFT_STRICT_HALF_SAVING.md#releasing-the-high-roots-with-a-determined-low-root)

Once the high roots are fixed, let
\(K_j=A_j\prod_{i\in I_j}x_i^{e_i}\).
On any remaining valid row the \(K_j\) are pairwise
coprime, and so are the \(z_j\). A letter with
\(T_t>2H_*^3\) has \(K_t>2H_*^2\). Modulo \(K_t\),
the linear relation
between the other two cofactors prescribes their
unit ratio, with the sign given by their positions
in the equation. If two pairs give this residue,
their cross determinant is divisible by \(K_t\)
and has absolute value at most \(2H_*^2\), so
their exact ratios agree. Reducedness recovers
both positive cofactors. The equation then
determines \(z_t\). Thus the remaining linear
fibre has size at most one. This also proves
the claim when there are no high roots.
The dimension descent has at most \(m+2\) stages.
Summing its covering exponents proves (MCB1).

## Bounded-group dictionaries with arbitrary fixed coefficients

**MCB2.** Use the root lists, raw cap \(A\), and
height-first conventions of (RSG1). Partition the roots
on each letter into groups of bounded size, each within
that letter. Give each group its joint dictionary cap
\(n^{w_G+o_X(1)}\), with \(w_G\) at most its summed raw
cap, and let \(w_T\) be the sum of these costs on \(T\).
Fix \(H_0>0\) and consider positive primitive rows
\(P+S=Q\) with \(Q\le X^{H_0}\). All three fixed
coefficients may be arbitrary within the complete
height cap. For every \(\eta>0\) there exist
\(C,L_0>0\) such that at any fixed \(0<L<L_0\),
if \(Q\ge n^C\), the count is
\[
 O\!\left(
       n^{(w-\max\{w_P,w_Q,w_S\})/2+\eta}
       X^\varepsilon
   \right)
                     \qquad(\varepsilon>0).
 \tag{MCB2}
\]
The thresholds can depend on \(A,H_0,\eta\)
and the group size bound. At fixed budget the
constant can also depend on the finite lists and
orders and \(\varepsilon\), uniformly in the
coefficients. This covers a fixed numerical
height floor \(Q\ge X^{h_0}\) for every fixed
\(h_0>0\) at sufficiently small \(L\).
There is no relative height gap or variable-height
floor on any prescribed letter.

Use the greedy \(n^{O(1)}\)-cofactor concentration
and the full-modulus release in (RSG1). For the
initial character moments on bounded groups, take
a lock \(F\) with \(\log_n F\to\infty\) and strip
it to the full divisor \(D\) of (RSG1). In the
\(2h_G\)-th moment choose the first \(h_G\) group
tuples at cost \(n^{h_G w_G+o_X(1)}\), with
\(h_G=\lceil A/b_G\rceil\) for the group's summed
raw cap \(b_G>0\). Collapse the remaining \(h_G\)
tuples into one product for each original position,
of height at most \(n^{2A}\). Apply the amplified
full-modulus release of (RSG1) to the collapsed
group monomial in its prescribed residue class.
There are boundedly many positions in the group,
and the fixed original position colouring gives
an exact monomial fibre of size one. Its count is
\(n^{o(1)}X^\varepsilon\), uniformly across groups,
since the root heights have one common bound and
the threshold in (GW) is independent of the fixed
exponents. Ordered decompositions of the collapsed
products are divisor-bounded at fixed moment
order. The norm therefore costs
\(w_G/2+o(1)/(2h_G)\) in base \(n\).
Hölder gives the opposite half bound, using
\(\sum_G(2h_G)^{-1}\le1\) and freezing zero-raw-cap
groups at subpower cost. Greedy packing now sets
aside boundedly many groups and hence boundedly
many original roots.

Call a variable product bounded if
its normalized logarithmic height is bounded
along the outer sequence. After passing to
subsequences, each whole variable product
is bounded or its normalized logarithmic
height tends to infinity. When applying
(MCB1) to a concentrated letter, the retained
list is bounded in number. Pass to a
subsequence on which each retained order is
bounded or tends to infinity. Move those of
bounded order into its \(n^{O(1)}\) cofactor.
The other orders tend to infinity and
(MCB1) counts the resulting terms at
\(n^{o(1)}X^\varepsilon\) cost. A fixed
complete letter or a recorded complement
can be absorbed into its coefficient.
A wholly bounded variable product is one
cofactor. These operations preserve the
pairwise coprimality of the variable
coordinates and the full positive
height floor on \(Q\).
[Radical-scale concentration and full-modulus release](ABC_RADICAL_SCALE_REAL_GAP.md#concentration-and-release-at-a-variable-gap)

If the thresholds did not exist, take a violating
sequence with \(L\to0\), normalized logarithmic
height of \(Q\) tending to infinity uniformly on
the counted rows, and limiting normalized count
exponent \(d\). A supported letter, meaning \(w_T>w-2d\),
has the concentrated description by the
first character-moment argument of (RSG1).
If \(w_T<d\) and another letter \(U\) is
supported, fix \(T\) completely. Every
lock on the remaining letter \(V\) of
normalized logarithmic size tending to
infinity releases \(U\) and must cost at
least \(d-w_T\). Thus \(V\) is concentrated
also. If the whole variable product of
\(T\) has normalized logarithmic size
tending to infinity, use it as a fixed
lock for the full-modulus release of both
other letters. If it is bounded, apply
(MCB1) to the two concentrated letters
with \(T\) fixed. Both give a contradiction.
It follows as in (RSG1) that
\(4d\le w+\min_T w_T\) and that there is
exactly one supported letter \(T_*\)
on the violating sequence.

Name the other two letters \(U,V\) with
\(w_V\le w_U\). Then \(w_V<d\) and
\(w_U+w_V<2d\). If both variable products
have normalized logarithmic size tending
to infinity, the two lock switches in
(RSG1) use only full-modulus releases and
give \(2d\le w_U+w_V\). If \(U\) is
bounded, fix \(V\) and apply (MCB1)
to \(T_*\) and the cofactor on \(U\).
If \(U\) has normalized logarithmic size
tending to infinity and \(V\) is bounded,
fix \(U\) when \(w_U<d\) and use the
same terminal. In the remaining case
fix \(V\) for the preliminary lock
argument. As in (RSG1), it gives a
bounded original group on \(U\) of cost at least
\(d-w_V\), whose powered product has
unbounded normalized logarithmic size.
Fixing its complement costs at most
\(w_U+w_V-d<d\). Release that bounded
group, \(T_*\), and the cofactor on \(V\)
by (MCB1). This is again a contradiction.
[Quarter inequality and the two lock switches](ABC_RADICAL_SCALE_REAL_GAP.md#counting-with-the-three-original-costs)

## Compressed fixed shifts without a unit shift modulus

**MCB3.** Fix \(m\ge0\) and \(B,A_0,H>0\). Use positive
coprime complete values \(U,V\) satisfying
\(\sigma_U U+\sigma_V V=t\ne0\), with fixed signs and
shift at each height and \(\max(U,V,|t|)\le X^H\).
One term has an integral representation
\[
            U=A u_0\prod_{i=1}^{m} x_i^{e_i},
       \qquad u_0,x_i\le n^B ,
\]
with \(A\) fixed and all \(e_i\) positive integers.
The other term has an anatomic root factorization
of total normalized raw cap at most \(A_0\) and
arbitrary fixed coefficient. Same-letter prime
overlaps in the compressed representation are
allowed. For every \(\eta>0\) there are fixed
\(C,L_0>0\) such that at any fixed \(0<L<L_0\),
if \(\max(U,V,|t|)\ge n^C\), the count is
\(O(n^\eta X^\varepsilon)\). In particular this
applies along any outer sequence with the maximum
normalized logarithmic height tending to infinity
and the compressed orders tending to infinity.
The constants are uniform in the fixed shift and
coefficients, and may depend on the original fixed
root lists at each budget.

Otherwise take a violating outer sequence with
\(L\to0\), \(\log_n\max(U,V,|t|)\to\infty\),
and fixed positive lower count exponent. Fix a
block on \(V\) with full product \(F\) and
\(\log_n F\to\infty\). The compressed variable
part of \(U\) lies in a prescribed unit residue.
Use the rebased (GW) release, with \(u_0\) as
the wide variable of height \(n^B\).
For fixed \(\ell\), an amplified exact monomial
\(\Phi\) is a positive integer of polynomial
height in \(X\). Each amplified root and the
amplified \(u_0\) divide \(\Phi\), regardless of
within-letter overlaps and regardless of whether
the orders exceed the modulus's effective
fullness parameter. Thus its exact fibre has
\(X^{o_X(1)}\) tuples by the divisor bound.
The individual coordinate products also admit
only \(X^{o_X(1)}\) ordered decompositions.
Amplification therefore gives release cost
\(n^{o(1)}X^\varepsilon\), as in (RSG1).

Every such lock on \(V\) must have positive
cost bounded below by the limiting count
exponent. Greedy packing leaves boundedly many
original roots and a cofactor of height
\(n^{O(1)}\) on \(V\). On the two bounded
root lists, pass to a subsequence on which each
order is bounded or tends to infinity.
Absorb the bounded-order factors into their
respective cofactors, increasing their
\(n\)-height bounds by only fixed constants.
The fixed shift is the third complete letter,
with its variable cofactor equal to one.
Signs can be arranged as a positive three-term
sum. Cross-letter coprimality follows from the
primitive pair and its shift. Since the maximum
complete height eventually exceeds \(2H_*^3\)
for the resulting fixed \(H_*=n^{O(1)}\),
(MCB1) counts this fibre at vanishing normalized
cost, a contradiction. This integer-cofactor
criterion strengthens the arbitrary-shift
branch of (RT2) in this discrete sector.
A correction of numerical height \(X^{o(1)}\)
alone need not have height \(n^{O(1)}\).
[Compressed fixed-shift contexts (RT2)](ABC_BALANCED_FIXED_SHIFT_REVERSE_TRANSFER.md#an-optimized-one-sided-cutoff-bound)

## Arbitrary-span progression occupancies

In the primitive quotient problem for a pair
\(X/2\le U<V\le X\) on a common unit progression
of step \(M\), write \(j=(V-U)/M\le n^s\) and
\(g=\gcd(U,V)\mid j\). At a fixed bin
\(g=n^{\chi+o(1)}\), the primitive equation
is
\[
           U/g+M(j/g)=V/g .
\]
The term \(j/g\) is one degree-one variable
of cap \(s-\chi\). If the original radical
caps are \(n^a,n^b\), let \(r_1,r_2\) be the
radical exponents of \(U/g,V/g\). Their
anatomic singleton dictionaries have costs
at most \(r_1,r_2\) up to arbitrary small
allowances. Since \(g\le n^s\), the complete
value \(V/g\) has a fixed positive numerical
height floor at sufficiently small budgets.
Apply (MCB2) with the prescribed coefficient
\(M\) on \(j/g\).

Restore \(g\) at normalized cost at most
\(y=\min(\chi,a-r_1,b-r_2)\) by counting its
radical outside the support of either quotient
and then its supported part. As usual one
can bin the normalized exponents on a fixed
fine mesh, with the mesh allowances absorbed
in \(\eta\). The three branches of the pair
bound give
\[
 \frac12\min\{r_1+r_2,r_1+s-\chi,r_2+s-\chi\}+y
       \le \frac12\min\{a+b,a+s,b+s\}.
\]
This extends the arbitrary-step close-pair
bound to arbitrary physical span, at arbitrary
valuation profiles. For a single unit block
with a common radical cap \(n^a\), its
occupancy exponent is at most
\(\min\{s,a/2,(a+s)/4\}+\eta\), by counting
pairs, taking a square root, and using the
trivial index-width bound.

For the occupancy of one fixed block, (MCB3)
gives a subpower bound. Fix \(s,a,H_0,h_0>0\).
Take integers in \([X^{h_0},X^{H_0}]\), of radical
at most \(n^a\), in one unit class modulo \(M\)
within a progression block of parameter width
\(n^s\). Choose the least eligible integer \(T\)
if the block is nonempty. For every other
eligible \(V\), write
\(V-T=M j\) with \(0<j\le n^s\) and
\(g=\gcd(V,T)\). Then \(g\mid T\), and the unit
class makes \(g\) coprime to \(M\), so \(g\mid j\).
Enumerate \(g\) at divisor cost. For each \(g\)
the primitive fixed-shift equation is
\[
                  V/g-M(j/g)=T/g.
\]
The variable part \(j/g\) is at most \(n^s\)
and can be used as the single integral cofactor
in (MCB3) with no high roots on that side.
The other term has an anatomic raw cap at most
\(a\) up to arbitrarily small allowances.
Since \(g\le n^s\), the normalized complete
height tends to infinity uniformly as \(L\to0\).
Thus for every \(\eta,\varepsilon>0\) the block
has \(O(n^\eta X^\varepsilon)\) eligible values
for all sufficiently small fixed budgets, with
no restriction on \(M\) beyond the unit class.
The quantitative pair bound above counts
across all classes and blocks, where a common
anchor \(T\) cannot be fixed at no cost.
[Anatomic factorization, Lemma 2.2](https://arxiv.org/html/2410.12234v2#S2)
· [Primitive quotient restoration](ABC_FIXED_SHIFT_STRICT_HALF_SAVING.md#progression-occupancies-with-highly-full-inputs)
