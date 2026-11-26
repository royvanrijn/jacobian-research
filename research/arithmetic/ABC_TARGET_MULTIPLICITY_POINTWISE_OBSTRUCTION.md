# Pointwise fibre bounds allow full target multiplicity on abstract additive labels

At every prescribed normalized multiplicity exponent between zero and
one half there is an abstract tripartite addition graph that realizes
the corresponding critical target-support exponent. Its labels can
be encoded by actual positive integers with distinct prime lists,
balanced heights, and a small joint radical. For every fixed label,
either free letter obeys the exponent in (FL1) after discounting
arbitrary polynomial-height prime pools. The addition relation is on
the abstract labels. The encoded integers serve as a counterexample
to deducing (FL8) from these pointwise and numerical constraints alone.
[Pointwise bounds and the target multiplicity estimate](ABC_FIXED_LETTER_RADICAL_CAPS.md)

## A graph with uniform intermediate multiplicity

Fix integers \(s\ge 2\) and \(0\le b\le s/2\), and set \(a=s-b\).
For each sufficiently large odd prime \(q\), let \(G\) be the space of
polynomials of degree less than \(a\) over \(\mathbb F_q\), and let
\(H\subseteq G\) consist of the polynomials of degree less than \(b\).
When \(b=0\) take \(H=\{0\}\). Use a copy of \(G\) for each of the three
parts. Put one edge on
\[
             (x,\ x+h,\ 2x+h),\qquad x\in G,\quad h\in H .
\]
The third label is the sum of the first two in \(G\), so any two
labels determine the third. There are \(q^{a+b}=q^s\) edges and all
\(q^a\) labels appear in each part. Fixing any one label leaves \(q^b\)
edges. More precisely, each of the other labels runs bijectively
through an affine coset of \(H\). For a fixed first label \(x\), this
follows by varying \(h\) in \(x+h\) and \(2x+h\). For a fixed second
label \(y\), the other labels are \(y-h\) and \(2y-h\). For a fixed
third label \(z\), they are \((z-h)/2\) and \((z+h)/2\). These formulas
also verify the assertion when \(b=0\).

## Actual prime-powered labels

Choose \(s\) distinct evaluation points \(t_1,\ldots,t_s\) in
\(\mathbb F_q\). Write
\[
      e_j=s+j-1,\qquad E=\sum_{j=1}^s e_j,\qquad
      L=\frac{3s}{E}=\frac{6}{3s-1}.
\]
Fix \(0<\vartheta<(\log 2)/(4E)\), and put
\(Q=D_s q\log q\) with a sufficiently large constant \(D_s\). For
every letter \(W\in\{A,B,C\}\), position \(j\), and symbol
\(u\in\mathbb F_q\), assign a distinct prime \(p_{W,j,u}\). Place
the \(C\)-primes in \([Q,(1+\vartheta)Q]\), and both sets of summand
primes in
\[
      [2^{-1/E}Q,\ (1+\vartheta)2^{-1/E}Q].
\]
The two displayed intervals are disjoint. The prime number theorem
gives enough primes for all the disjoint lists when \(q\) is large
at fixed \(s\) and \(\vartheta\), since each interval then contains
asymptotically a positive constant times \(D_s q\) primes. The
encoding on part \(W\) is
\[
               W(f)=\prod_{j=1}^{s}p_{W,j,f(t_j)}^{\,e_j}.
\]
The evaluations determine \(f\), since two polynomials of degree
less than \(a\le s\) agreeing at \(s\) distinct points are equal.
Hence each encoding and each radical encoding is injective on its
own part. All three integers on an edge are pairwise coprime.

Put \(X=(1+\vartheta)^E Q^E\). Every encoded \(C\) lies in
\([Q^E,X]\), and every encoded \(A,B\) lies in \([Q^E/2,X/2]\).
In particular \(X/2<C\le X\), and the separate height ranges
for \(A+B\) and \(C\) both lie in \([Q^E,X]\). Moreover,
\[
 \operatorname{rad}(ABC)
   \le 2^{-2s/E}(1+\vartheta)^{3s} Q^{3s}
   < Q^{3s}
   \le C^L .
\]
At each fixed \(s\) we have
\(\log\operatorname{rad}(W)/\log C=L(1/3+o_q(1))\)
uniformly in the labels. The individual prime masses in base
\(X^L\) equal \(1/(3s)+o_q(1)\). Consequently, after taking
\(q\to\infty\) at each fixed \(s\) and then \(s\to\infty\),
each letter has the exact-valuation profile
\((1/3)\operatorname{Unif}[2,4]\) and numerical height one.
The notation \(\operatorname{Unif}[2,4]\) denotes the uniform
probability measure on that interval. The height limit precedes
the budget limit.

## Uniform bounds for arbitrary prescribed prime pools

Fix one encoded complete label \(T\) and prescribe a pool of
primes dividing a positive integer \(P\le X^B\), for fixed
\(B>0\). The pool may depend on \(T\) and on \(q\). Let \(W\) be
either free letter. Every prime appearing in its labels is at
least \(Q_-=2^{-1/E}Q\). Thus the pool contains at most
\[
                      \frac{B\log X}{\log Q_-}=O_{s,B}(1)
\]
of the primes used by the encoding. The number of possible subsets
of captured primes is bounded in terms of \(s,B\) alone. Each subset
with \(k\) captured primes fixes \(k\) evaluation coordinates of the
free label, provided the subset occurs. Two primes from the same
position cannot occur in one such label.

The projection from \(H\) onto any \(k\) distinct evaluation
positions has rank \(\min(k,b)\), by the Vandermonde matrix or
polynomial interpolation. Its translate on a fixed-label fibre
has the same rank. Thus any one record of \(k\) captured primes
leaves at most \(q^{(b-k)_+}\) rows. Since \(b\le s/2\) and
\(k\le s\),
\[
                      (b-k)_+\le (s-k)/2 .
\]
If the outside radical on \(W\) is at most \(X^{L\ell_W}\),
its \(s-k\) remaining primes give \(Q_-^{s-k}\le X^{L\ell_W}\).
For large \(q\), \(q\le Q_-\), so the number of rows for this record
is at most \(X^{L\ell_W/2}\). Sum over the boundedly many captured
subsets. By applying this argument separately to two possibly
different pools on the two free letters, the common restriction
has at most
\[
             O_{s,B}\!\left(X^{L\min(\ell_U,\ell_V)/2}\right)
\]
rows. This is the pointwise exponent of (FL1), uniformly in the
fixed label and in the pools. It also gives the bound after
recording full powers on the pools, since the order of each
prime used in the encoding is prescribed.

## Consequence at each multiplicity level

The maximum \(C\)-radical is at most
\(Z=(1+\vartheta)^s Q^s\). One can equivalently use \(Z=Y^{1/M}\)
with \(Y=Z^M\) and any fixed \(M\ge\max_j e_j\), since the
least target repair sends the encoded \(C\) to
\(\operatorname{rad}(C)^M\). This use of the repair changes the
target cutoff only. It makes no claim about primitive quotients
of an actual sum.

The \(q^a\) distinct target radicals each carry \(q^b\) rows.
Since \(\log Z=(s+o_q(1))\log q\), the exponents of the target
support, multiplicity, and total count are respectively
\[
                    1-b/s,\qquad b/s,\qquad 1.
\]
For every real \(u\in[0,1/2]\), choose \(b/s\to u\) as
\(s\to\infty\). In particular this construction includes every
interior range of multiplicities below the square-root fibre bound
as well as profiles approaching that bound. For a fixed \(b>0\)
and a threshold \(v<b/s\), eventually every support has
multiplicity at least \(Z^v\). Given any proposed saving
\(\sigma>0\), one may take \(v\) within \(\sigma/2\) below \(b/s\).
Then \(1-b/s>1-v-\sigma\), contradicting the analogue of (FL8) on
these abstract edges. At \(b=0\) the same assertion holds at
threshold \(v=0\).

## Independence from a prescribed finite congruence test

The labels can be assigned with at most one true integer sum on
the abstract edges. Keep the summand encodings fixed. For the
target, independently permute the assignments of the \(q\)
available primes to the \(q\) symbols at each position. For a given
edge, the value \(A+B\) is fixed. Unique factorization and the
disjoint prime lists give at most one choice of target primes
that realizes this value. The independent permutations assign
that tuple to the edge's target label with probability at most
\(q^{-s}\). Thus the expected number of true sums is at most one.
Some assignment realizes this bound, and all the pointwise
arguments hold for every assignment.

The same observation works after imposing the integer congruence
\(A+B\equiv C\pmod J\) on every edge for any fixed \(J\ge2\).
Choose all list primes congruent to \(1\bmod J\), which is possible
by the prime number theorem for arithmetic progressions at fixed
\(J\). Use the original summand centres \(2^{-1/E}Q\), use the
same centre for the target lists, and multiply the encoded target
by \(2\). The target then lies in \([Q^E,X]\) as before, and the
summands and target have the asserted congruence. The factor \(2\)
is coprime to all primes in the lists. Replace \(L\) by
\(L'=(3s+1)/E\); then the combined radical is
\(O_{s,J}(Q^{3s})<C^{L'}\) for large \(q\). The target cutoff
changes by at most a factor depending on \(s,J\) at each fixed
budget, and \(\log\operatorname{rad}(W)/(L'\log C)\to s/(3s+1)\)
as \(q\to\infty\). The half-cap proof still applies with \(L'\), counting just the
list primes when forming the set of captured positions.
These changes preserve the limiting balanced profile and the
multiplicity exponents. Randomly permuting the target lists
as above preserves the congruence while making at most one edge
an actual sum. A fixed finite list of congruence tests can be
imposed by taking \(J\) to be their least common multiple.

The earlier support-labelled model in (DRR5)–(DRR7) and its
following obstruction distinguish local from joint divisor
frequencies. Here the quantities specified are the support size
and the multiplicity on every complete-label fibre, including the
uniform pool bounds at intermediate multiplicities. The model
isolates the arithmetic input that is still needed in applying
(FL8) to actual primitive sums.
[Reciprocal records and the support-labelled obstruction](ABC_DIVISOR_RECORD_RECIPROCITY.md#reciprocal-records-on-a-large-fixed-value-fibre)
