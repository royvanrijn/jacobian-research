# Release with two compressed letters

If two letters are each a fixed coefficient times boundedly
many powered roots and a small integral cofactor, a primitive
sum gives a subpower count even when the third letter has an
arbitrary valuation profile. A fixed factor on either
compressed letter need not be full. The full modulus used in
the proof comes from the original roots of the third letter.

## Two compressed terms

Fix \(A,B,H,\theta>0\) and a nonnegative integer \(m\).
Put \(n=X^L\), and fix the finite root lists and their
orders before taking \(X\to\infty\) at each budget \(L\).
Count positive primitive rows \(A_1+A_2=A_3\) at height
\(Q=\max_j A_j\le X^H\). On two prescribed distinct
roles \(U,V\), suppose the original complete values have
representations
\[
 U=F_U u_0\prod_{i\in I_U}x_i^{e_i},
 \qquad
 V=F_V v_0\prod_{i\in I_V}x_i^{e_i},
 \qquad |I_U|+|I_V|\le m ,
\]
where the positive \(F_U,F_V\) are fixed, the orders are
positive integers, and \(x_i\le n^B\). The cofactors
\(u_0,v_0\) are positive integers and can share primes
with the roots on their own letter. On the third role
\(W\), use an anatomic factorization into a fixed
positive coefficient times powered original roots.
The sum of the normalized raw root caps on \(W\)
is at most \(A\).

**DTC1.** For every \(\eta>0\) there are \(C,L_0>0\)
such that for each fixed \(0<L<L_0\), if
\[
          u_0,v_0\le n^B,\qquad Q\ge n^C ,
\]
the count is \(O(n^\eta X^\varepsilon)\) for every
\(\varepsilon>0\). There are also \(\xi,L_1>0\)
such that, at each fixed \(0<L<L_1\), the same count
holds under
\[
       u_0,v_0\le n^{B L^{-\xi}},
       \qquad Q\ge n^{L^{-\theta}} .
                                                               \tag{DTC1}
\]
The thresholds depend on the fixed cap bounds and
\(\eta\). The height-limit constants may also depend
on the finite lists, fixed budget and \(\varepsilon\).
They are uniform in the prescribed coefficients.
A fixed positive numerical height floor on \(Q\)
suffices for the second assertion by taking
\(0<\theta<1\).

### Releasing the compressed ratio

At each fixed budget refine the original roots on
\(W\) into dyadic height boxes at subpower cost.
Use prime-role colouring at the same cost, so a prime
in a variable numerator on one role never occurs in
a variable denominator on the other role, even across
different retained rows. This follows by random
colouring of the prime support, since a row uses
\(o_X(\log X)\) primes at polynomial height. The
fixed coefficient primes can be assigned their own
roles beforehand. They cannot conflict on a primitive
row containing the fixed factors.

Fix a block of original powered roots on \(W\) with
full product \(F\), and write \(w=\log_n F\).
The signed equation makes the two other complete
values have a fixed unit ratio modulo \(F\).
Discard the positions in this block whose orders are
less than \(\lfloor w/(2(A+1))\rfloor\).
They contribute at most half its normalized
logarithmic height, up to the dyadic allowances.
The remaining divisor \(D\mid F\) is \(k\)-full with
\(k=\lfloor w/(2(A+1))\rfloor\) and
\(\log_n D\ge w/2+o_X(1)\).
Original roots have fixed orders at the height limit,
so this trimming can use the limiting height of the
chosen block at that budget.

Absorb \(F_U,F_V\) into the unit target. The signed
monomial to which (GW) is applied is the ratio of
the two variable parts, with one rational cofactor
\(u_0/v_0\). Its narrow roots have height \(D^{O(1/k)}\)
when \(w\) tends to infinity. In the fixed-cofactor
case the wide exponent is \(O(1/k)\). In the
stretched-cofactor case require \(w\ge L^{-2\xi}\).
Then the wide exponent is
\(O(L^{-\xi}/w)=O(k^{-1/2})\), which tends to zero.
The estimates remain valid with the constants
multiplied by any fixed amplification order.

Take coordinatewise products of \(\ell\) tuples in
one fixed record. Their signed wide cofactor has no
cross-role numerator-denominator cancellation by
the colouring. Every amplified narrow root and
each amplified cofactor divides the corresponding
numerator or denominator of the exact monomial
value. All exact values have polynomial height
in \(X\) for fixed \(\ell\), because they are products
of ratios of the original complete terms divided
by their fixed coefficients. By the divisor bound,
the exact fibre and the ordered decompositions of
the amplified coordinate products are subpower.
This argument allows overlaps between a narrow
root and the cofactor on the same letter, and does
not require its order to exceed \(k\).

Apply (GW) to the amplified tuples with target
raised to the \(\ell\)-th power. Its main exponent
in base \(n\) is \(O_m(\log_n D/k)=O_{m,A}(1)\).
Taking \(\ell\)-th roots and choosing an arbitrarily
large fixed \(\ell\) gives release at vanishing
normalized cost in the outer limit. The exact
ratio of the complete free terms recovers both by
primitivity, and hence the third through the original
signed equation. All modulus data and amplification
orders are fixed for the height limit.
[Full-modulus count (GW)](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-bounded-number-of-small-roots)
· [Release at a varying normalized scale (RSG1)](ABC_RADICAL_SCALE_REAL_GAP.md#concentration-and-release-at-a-variable-gap)

### Concentrating the third letter

For the fixed-cofactor assertion, if the stated
thresholds did not exist take a sequence of
violations with count exponent bounded below by
\(\eta>0\), \(L\to0\), and \(\log_n Q\to\infty\).
The ratio release implies that every block on \(W\)
with normalized logarithmic height tending to
infinity would cost at least \(\eta/2\) to fix at
late budgets. Set aside the original roots of raw
cost at least \(\eta/8\). Pack the remainder greedily
into blocks of cost at most \(\eta/4\).
The list set aside and the number of blocks are
both bounded in terms of \(A,\eta\). Each cheap
block has bounded normalized logarithmic height
along the sequence, since an unbounded block would
give a violating subsequence. Thus the remaining
part of \(W\) is one integral cofactor of height
\(n^{O(1)}\). Among the boundedly many retained
roots on all roles, pass to a subsequence with
orders bounded or tending to infinity. Absorb the
bounded-order roots into their own cofactors.
The other orders tend to infinity, and
\(Q>2H_*^3\) for the resulting bound \(H_*=n^{O(1)}\)
on all roots and cofactors. The multiple-cofactor
boundary bound (MCB1) counts the original rows
at vanishing normalized cost, a contradiction.

For the stretched-cofactor assertion, choose
\(\xi>0\) sufficiently small in terms of
\(A,B,m,\eta,\theta\). Along a violating sequence
of count exponent at least \(\eta\), the ratio
release forbids any block of cost at most
\(\eta/4\) on \(W\) whose log-height is at least
\(L^{-2\xi}\) in base \(n\). The same packing now
leaves a cofactor on \(W\) of height
\(n^{O(L^{-2\xi})}\) and a bounded list of original
roots. The initial lists on \(U,V\) are also
bounded. Apply (TW1) to all three integral
cofactors and these original roots, taking
\(2\xi\) strictly below its width exponent for
the worst possible retained number of roots.
This absorbs the fixed constants in the
cofactor heights. The lower height threshold
\(Q\ge n^{L^{-\theta}}\) is precisely the
condition for (TW1). The resulting
subpower bound is again a contradiction.
There are no fixed numerical-height or
effective-fullness assumptions on \(F_U,F_V\).
[Boundary release with multiple small cofactors (MCB1)](ABC_MULTIPLE_SHORT_COFACTOR_BOUNDARY.md#releasing-the-high-roots-and-all-three-cofactors)
· [Finite order gaps (TW1)](ABC_THIN_WINDOW_ARBITRARY_SHIFT_RELEASE.md#a-gap-among-finitely-many-orders)

## A joint two-divisor record across the balanced range

**DTC2.** Consider actual primitive triples of
height \(X/2<C\le X\) with bounded individual
radical caps in base \(n\). Suppose on two
prescribed roles \(U,V\), each row admits
divisors \(F_U\mid U,\ F_V\mid V\) from a
joint dictionary of size at most \(n^f\), with
\[
               U/F_U,\ V/F_V\le n^{L^{-\xi}}
\]
for the sufficiently small fixed positive
\(\xi\) of (DTC1). Then the row count is
\(O(n^{f+\eta}X^\varepsilon)\). In particular
this excludes balanced saturation on any
sector admitting such a dictionary at
\(f<1/3\) by a fixed margin. The assertion
uses no lower bound on complete-letter
multiplicity and applies throughout the
range of \(\alpha\) in (RB3), including
occurrence.

Fix the two divisors and apply (DTC1) with
no powered roots on either compressed
side. Anatomic factorization of the third
role fixes any exceptional coefficient and
the finite root-height boxes at arbitrarily
small normalized cost. The original total
maximum supplies the (DTC1) height floor.
Sum over the joint catalogue and decrease
the internal allowances before choosing
\(\xi\). The counted rows can retain all
original radical and least-repair
restrictions. There is no introduction of
a common multiplier.
[Balanced dictionary exponents (RB3)](ABC_BALANCED_DICTIONARY_RESIDUE_CELLS.md#equality-of-the-three-complete-dictionary-exponents)
· [Anatomic factorization, Lemma 2.2](https://arxiv.org/html/2410.12234v2#S2)

The same proof allows a fixed number of
thin scaled-degree windows remaining on
each of the two roles once the two records
are fixed. In each window put the product
of its original roots to the minimum
original order. If the window widths are
\(O(L^{1-\xi'})\) for some \(\xi'<\xi\),
the exponent differences put at most
\(n^{O(L^{-\xi'})}\) in each combined
integral correction. The grouped roots
have height \(n^{O(1)}\), so (DTC1)
applies with the corresponding fixed
root number and cap bounds. The catalogue
must charge the entire outside record on
the two roles together. Taking two
complete letters as that record costs
the whole row exponent, since two
complete values determine the row.
The saving requires a cheaper joint
record. The result therefore
does not extract an inexpensive
two-divisor record from arbitrary
balanced spectra.
