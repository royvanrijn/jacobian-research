# Real-gap release at the radical scale

For singleton root dictionaries the real height gap in (SHF4)
can shrink with the budget. The proof first uses unit
character moments under a modulus of size \(n^{O(1)}\)
to put the variable part of a concentrated letter into
boundedly many root coordinates and a cofactor of height
\(n^{O(1)}\). Both the real and full-modulus release
estimates can then be applied with rescaled height parameters.
[Grouped short-letter estimate (SHF4)](ABC_FIXED_SHIFT_STRICT_HALF_SAVING.md#grouped-dictionaries-with-a-large-coefficient-on-the-short-letter)

## Statement

Fix \(A,H,h_0,\zeta>0\) with \(0<\zeta<h_0\le H\).
At each fixed budget \(L>0\), put \(n=X^L\) and take
the height limit with the finite root lists and their
positive integer orders fixed. Consider actual positive
primitive rows \(P+S=Q\). Write each letter as a fixed
positive integer coefficient times a product of powered
original roots. The roots on all original positions
of a row are pairwise coprime. The two tall letters
\(P,Q\) have heights in \([X^{h_0},X^H]\), and their
fixed coefficients are at most \(X^{h_0-\zeta}\).
The coefficient of \(S\) may have arbitrary height
within the complete height bound \(S\le X^H\).

For each original root \(z_i\), let \(n^{b_i}\) be its raw
cap and \(n^{w_i+o_X(1)}\) a cap on its dictionary
cardinality, where \(0\le w_i\le b_i\) and \(\sum_i b_i\le A\).
Write \(w_T\) for the sum of the singleton dictionary costs
on letter \(T\), and \(w=w_P+w_Q+w_S\).

**RSG1.** For every \(\eta>0\) there are \(C,L_0>0\) such
that, for fixed \(0<L<L_0\), if
\[
                         S\le Q n^{-C},
\]
the row count is
\[
       O\!\left(
          n^{(w-\max\{w_P,w_Q,w_S\})/2+\eta}X^\varepsilon
          \right)
                   \qquad(\varepsilon>0).
 \tag{RSG1}
\]
The thresholds depend on the fixed cap and height margins and
\(\eta\). The constants in the height limit may also depend on
the chosen budget, finite degree lists and \(\varepsilon\).
They are uniform in the fixed coefficients within their stated
height ranges. The conclusion is for singleton dictionaries.

## Concentration and release at a variable gap

At fixed \(L\), refine all root sizes dyadically at
\(X^{o_X(1)}\) cost so that \(\log_n z_i=b_i+o_X(1)\)
within their boxes, decreasing the stated raw and dictionary
caps as appropriate. Work on height subsequences to give
these exponents fixed limits at each budget. Positions of
raw cap zero are fixed at height-subpower cost. The orders
are fixed in this step, so the normalized logarithmic size
of any selected product of powered roots differs from
\(\sum_i e_i b_i\) for those positions by \(o_X(1)\).
Pass also to a common colouring of all variable primes by
their original positions, at \(X^{o_X(1)}\) loss. Indeed each
row's probability of satisfying independent uniform colours
is \(J^{-\omega(PQS)}=X^{-o_X(1)}\), where \(J\) is the
number of colours at this fixed budget. A prime of a
variable root now stays in its position throughout every
product argument.

First fix any block of variable positions on one letter
whose full product \(F\) is larger than \(n^{2A+1}\).
The other two letters have a prescribed unit ratio modulo
\(F\), since the original triples are primitive. The cost
for this congruence is at most half the sum of the original
dictionary costs on these two letters. To prove this, apply
unit character orthogonality to the signed variable
monomial, absorbing the fixed coefficients into the
target. For position \(i\), use the \(2h_i\)-th moment
with \(h_i=\lceil A/b_i\rceil\). Products of \(h_i\) roots
are at most \(n^{2A}\), and the number of \(e_i\)-th roots
of unity modulo \(F\) is at most \((2e_i)^{\omega(F)}\).
For an odd prime power the unit group is cyclic, while
for a power of two it is the product of a cyclic group
with a group of order at most two. This proves the
bound, which is \(X^{o_X(1)}\) uniformly for \(F\le X^H\)
and fixed \(e_i\). Given one collapsed product, every
possible residue for the second product determines at
most one integer in its range. Ordered product
factorizations cost \(X^{o_X(1)}\) by the divisor bound.
Thus the \(2h_i\)-norm is at most
\(n^{w_i/2+o_X(1)}\). Since
\(\sum_i (2h_i)^{-1}\le1\), Hölder gives the asserted
half cost. The argument applies also to a single free
letter in a prescribed unit residue, and applies when
some original positions have already been fixed.
[Character moments for an arbitrary unit modulus (RO3)](ABC_RESIDUE_OCCURRENCE_ROOT_PROJECTION_OBSTRUCTION.md#a-residue-bound-at-several-depths)

We will use the following greedy consequence. Suppose
that along a sequence of boxes every block with
unbounded normalized logarithmic size, tending to
infinity on a subsequence, is forbidden to have cost
below some fixed \(\beta>0\). Set aside the singleton
positions costing at least \(\beta/4\). Greedily pack
the others into blocks of cost at most \(\beta/2\).
All completed blocks except possibly the last cost
more than \(\beta/4\), so their number is bounded in
terms of \(A,\beta\). Each of these cheap blocks must
have bounded normalized logarithmic size along the
sequence. Otherwise one of them would tend to infinity
on a subsequence and violate the assumed cost floor.
Thus outside boundedly many original roots the entire
variable product is \(n^{O(1)}\). If the whole variable
product has normalized logarithmic size tending to
infinity, at least one of the set-aside singletons
does as well after a subsequence.

Next consider a fixed block on the omitted letter with
full product \(F\) and \(u=\log_n F\to\infty\) in the
outer limit. One can take \(u\) to be the limiting
exponent in its fixed-budget box. Let
\(k=\lfloor u/(2(A+1))\rfloor\). Drop the positions with
orders below \(k\), losing at most \(kA+o_X(1)\) from the
normalized logarithmic size. The remaining divisor \(D\)
is \(k\)-full and has \(\log_n D\ge u/2+o_X(1)\).
Hence boundedly many roots of height at most \(n^{O(1)}\)
and a cofactor of height at most \(n^{O(1)}\) meet the
narrow and wide parameters of (GW) based on \(D\).
In its notation the narrow height cap is \(D^{C_1/k}\)
and the wide height exponent can be taken \(C_2/k\),
with fixed constants \(C_1,C_2\). Here \(k\to\infty\),
and \(\log_n D/k\) stays bounded above.

The exact monomial for release is the ratio of the
variable parts of the two free letters, or a single
variable part if the other free value is fixed. Move
already recorded factors and the original fixed
coefficients into its target. Use the selected
original roots with their original signed orders and
use the signed product of the remaining factors as
the one rational cofactor. The fixed position
colouring separates the prime support of every
retained root from that of the cofactor and all other
retained roots. Thus the prime valuations of the
exact monomial recover the retained roots regardless
of the sizes of their orders relative to \(k\), and
then recover the exact cofactor. A fixed exact ratio
of the complete free values determines both by
primitivity, and the third letter is then determined
by the equation. If one free value was fixed, the
exact monomial gives the other value directly.

For a fixed amplification order \(\ell\), take
coordinatewise products of \(\ell\) such tuples
in a fixed record. Every cofactor prime belongs
throughout to the same signed side, so there is no
cancellation between numerator and denominator in
these products. The same position colouring gives
exact monomial fibre size one for the amplified
tuples. The preimage count of each amplified tuple
is \(X^{o_X(1)}\) by factoring its coordinate
products, including the separate numerator and
denominator products. Apply (GW) to these tuples
with target raised to the \(\ell\)-th power and
height constants multiplied by \(\ell\).
This bounds their normalized count exponent by
\(O(1)\), where the main constant comes from
\((m^2+m)\log_n D/k\) and is independent of
\(\ell\). Taking \(\ell\)-th roots and choosing
arbitrarily large fixed \(\ell\) gives a release
cost tending to zero in base \(n\).
This argument does not use the overlap-compatible
fibre bound that assumes orders comparable to \(k\).
[Bounded-root full-modulus count (GW)](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-bounded-number-of-small-roots)

For real release suppose \(P/Q=1+O(n^{-c_L})\) with
\(c_L\to\infty\). On fixed records where the relevant
free parts consist of boundedly many original roots
of height \(n^{O(1)}\) and one cofactor of height
\(n^{O(1)}\), change the base to \(Y=n^{c_L}\).
The budget parameter in (RUW1) becomes \(1/c_L\),
and the real precision parameter is \(1\). The
cofactor has \(Y\)-height exponent \(O(1/c_L)\),
so it lies inside any fixed positive wider-height
allowance for small budgets. The position colouring
and fixed-record ratio argument give exact fibre
size one. Applying (RUW1) with arbitrarily small
fixed chart allowance gives release at
\(n^{o(1)}X^\varepsilon\) cost in the outer limit.
The new parameters and all orders are fixed before
the height limit at each budget.
[Unequal-height real-strip release (RUW1)](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-real-strip-with-unequal-root-heights)

## Counting with the three original costs

Suppose (RSG1) fails. Diagonalize failures with
\(L\to0\), a fixed positive excess over the claimed
coefficient, and \(S/Q\le n^{-c_L}\) with
\(c_L\to\infty\). Pass to subsequences fixing limits
of the normalized costs and the row exponent \(d\).
Use the same \(w_T,w\) for their limits. In particular
the variable products on both tall letters have
normalized logarithmic size tending to infinity by
the fixed coefficient margins. Call a letter
supported if \(w_T>w-2d\). Fixing any block on a
supported letter with normalized logarithmic size
tending to infinity allows the first character
moment bound above. Consequently such a block must
cost at least \(d-(w-w_T)/2>0\). The greedy
consequence gives the bounded-root and
\(n^{O(1)}\)-cofactor description on every supported
letter.

If a letter \(T\) has \(w_T<d\), no other letter can
be supported. Indeed, suppose \(U\ne T\) is
supported and \(V\) is the remaining letter. Fix
\(T\) completely at cost at most \(w_T\). Any lock
on \(V\) whose normalized logarithmic size tends
to infinity releases \(U\) by the full-modulus
argument. Therefore such locks on \(V\) cost at
least \(d-w_T\). The greedy consequence also
concentrates \(V\). If \(T\) is tall, its whole
variable product is now a fixed lock of unbounded
normalized logarithmic size, allowing both \(U,V\)
to be released at vanishing additional cost. If
\(T\) is short, \(U,V\) are tall and the real
release counts the two together. Both cases
contradict \(w_T<d\).

It follows that \(4d\le w+\min_T w_T\). Otherwise
a cheapest letter costs less than \(d\), and the
two other costs sum to more than \(2(w-2d)\),
so at least one other letter is supported. In
particular \(d\le w/3\). A failure of the claimed
bound makes a most expensive letter \(T_*\)
supported. Three supported letters would imply
\(d>w/3\). Two supported letters would leave a
third with cost below \(4d-w\le d\), which was
excluded. Thus \(T_*\) is the only supported
letter.

Name the other two letters \(U,V\) with \(w_V\le
w_U\). Then \(w_U+w_V<2d\) and \(w_V<d\).
First suppose both variable products on \(U,V\)
have unbounded normalized logarithmic size
tending to infinity. Fix \(V\) completely. Any
such lock on \(U\) releases \(T_*\), so costs
at least \(d-w_V\). The greedy consequence
finds an original singleton on \(U\) of cost
\(\beta_U\ge d-w_V\) whose powered factor has
normalized logarithmic size tending to infinity.
Its complement has cost \(w_U-\beta_U<d\).
After fixing that complement, a lock on \(V\)
of normalized logarithmic size tending to infinity
releases the selected root on \(U\) and the
concentrated free part of \(T_*\) together. Thus
such locks on \(V\) have cost at least
\(d-w_U+\beta_U>0\). There is a singleton on
\(V\) with unbounded normalized logarithmic size
as well. Write its cost as \(\beta_V\).
Fixing the selected singleton on \(U\) and the
complement of the selected singleton on \(V\),
then using that \(U\) lock to release the other
variable parts, gives \(d\le \beta_U+w_V-\beta_V\).
Switching gives \(d\le\beta_V+w_U-\beta_U\).
Their sum contradicts \(w_U+w_V<2d\).

Otherwise, after passing to a subsequence, one
of \(U,V\) has bounded normalized variable
logarithmic height. It must be \(S\), since
both tall variable products have unbounded size.
So \(T_*\) is tall. If the other tall letter
costs less than \(d\), fix it completely and
use the real release on \(T_*\). Otherwise that
other tall letter is \(U\), since \(w_V<d\).
The first locking step with \(V\) fixed still
finds the singleton on \(U\) with complement
cost \(w_U-\beta_U<d\). Fix the complement and
use the real release for the ratio of the
two tall free parts. Both alternatives
contradict \(d\). This proves (RSG1), since a
failure for every \(C,L_0\) at fixed allowance
gives exactly the diagonal sequence used above.

## Consequence for arbitrary-step pair counts

In the primitive quotient reduction for close
pairs on a unit progression, write
\(g=\gcd(U,V)\), \((V-U)/M=gJ\).
Then \(S=MJ\) has the same height ratio to \(Q=V/g\)
as \(V-U\) has to \(V\). Its fixed coefficient
is \(M\) and its single root \(J\) has the
usual index-width cap. Use anatomic singleton
factorizations on the tall primitive quotients
with coefficient heights tending to zero in base
\(X\). These heights leave a fixed margin to
the tall complete heights. Hence the
arbitrary-step close-pair argument using (SHF4)
also works for span at most \(Xn^{-C}\), with
\(C\) depending on the required allowance and
fixed raw caps. The restoration of \(g\) and
all radical charges are unchanged.
[Anatomic factorization, Lemma 2.2](https://arxiv.org/html/2410.12234v2#S2)
