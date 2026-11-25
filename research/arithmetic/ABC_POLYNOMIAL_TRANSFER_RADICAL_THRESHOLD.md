# Radical thresholds for polynomial transfers

Use an affordable fixed record as in
[the projective packing theorem](ABC_MONOMIAL_TRANSFER_HEIGHT_AND_FERMAT_GEOMETRY.md#a-uniform-small-output-bound-for-rational-charts).
Write \(n=X^L\), with \(N=n^{r+o_X(1)}\) actual primitive rows at
height at most \(X\) with distinct targets. Assume the uniform input
radical and affordability margins of that theorem and that \(r\)
stays in a fixed compact subinterval of \((0,\infty)\).
Fix \(m\ge2\). In slot \(i\) choose distinct original letters
\(W_i,Z_i\), put \(x_i=Z_i/W_i\), and take a rational-coefficient
polynomial \(P(\boldsymbol x)\) of positive degree in every slot.
Degrees and polynomial coefficient-height bounds may be fixed
separately at the chosen \(L\) before the height limit. For
\(P(\boldsymbol x)=p/q\) in lowest terms with \(q>0\), call
the tuple eligible if \(p(p-q)\ne0\) and
\[
                     \operatorname{rad}(|pq(p-q)|)\le n^b.
\]
For every fixed \(\Delta>0\) and all sufficiently small fixed
\(L\), if \(0\le b\le4r-\Delta\), there are
\(O(n^{mr-\kappa+o_X(1)})\) eligible tuples for some fixed
\(\kappa>0\). The small-budget and density thresholds depend on
\(m,\Delta\) and the uniform fibre parameters. The inner losses
can depend on the polynomial degree and height bounds at fixed
\(L\). There is no output-height restriction.

## Counting with a pure denominator

Regularize the original tag multiplicities in the eligible graph
by dyadic classes, at subpower cost. In slot \(i\) write
\(n^{q_i+o_X(1)}\) for the tag degree cap, \(a_i=r-q_i\) for
the tag dictionary exponent cap, and \(s=\sum_i a_i\). The
lower and upper multiplicity bounds differ by only a constant
factor before restriction to the eligible graph. Bin the
pure output radical exponent near \(c\), in base \(n\).
An eligible graph of exponent arbitrarily close to \(mr\)
requires, up to arbitrarily small fixed allowances,
\[
 \begin{split}
 r-a_i&\le (b-c)/4,\\
 r+s-a_i&\le
        \min\{(b+c)/4,\ (b-c)/2,\ b/3\}\qquad(1\le i\le m).
 \end{split}                                                       \tag{FR1}
\]
Adding the two inequalities with right sides \((b-c)/4\) and
\((b+c)/4\) and averaging over \(i\) yields
\[
     b\ge4r+2(1-2/m)s,\qquad b\ge3(r+s-\min_i a_i).                  \tag{FR2}
\]
The second bound follows from the \(b/3\) upper exponent.

To prove (FR1), fix \(K\) large in terms of the required output
counting allowance. Discard the outputs of height at most \(n^K\)
at power-small cost by the tuple-chart theorem with any fixed
\(0<\rho<r/2\). Each ratio of distinct original letters is a
Möbius function of the primitive cofactor direction, with
polynomial coefficient height. For each slot discard also the
tuples where the specialized polynomial loses its positive degree
in that slot. The exceptional complementary tuples number
\(O_P(N^{m-2})\) by the Cartesian zero-set bound.
[The tuple-chart bound](ABC_MONOMIAL_TRANSFER_HEIGHT_AND_FERMAT_GEOMETRY.md#tuple-charts-and-the-pure-term-pattern)

The remaining heights are between \(n^K\) and \(X^{C(L)}\).
Bin the output height exponent in base \(n\) in intervals
\([h,h+\xi]\), \(h\ge K\), and bin \(\log_n\operatorname{rad}(q)\)
in \([c,c+\xi]\). The number of classes is fixed at the chosen
\(L\). Refine the height class into ordinary dyadic shells with
upper height \(Y\), and use the fixed relative budget
\(\lambda'=b/h\le b/K\) in that class. One may first round
\(b\) upwards on a grid of mesh \(\xi\). The budget on the
shell is at most \(Y^{\lambda'}\le n^{b+b\xi/K+o_X(1)}\).
Take the small fixed \(\xi\) and sufficiently large \(K\)
after choosing the desired output allowance. An error
\(\eta'\lambda'\) in the relative exponent contributes at
most \(\eta'(b+b\xi/K+o_X(1))\) in base \(n\).
The arbitrary height exponent allowances can also be taken
arbitrarily small in base \(n\), because the height exponents
in that base are bounded at fixed \(L\). Summing the dyadic
shells costs an inner subpower factor. Each signed output
gives a positive primitive triple after permutation, with
at most six signed ratio outputs per reordered triple.

Fix the tag \(W_i\) and all complementary complete rows. If
\(d_i\) is the degree in slot \(i\) and \(J\) clears coefficient
denominators, the unreduced denominator \(J\prod_j W_j^{d_j}\)
is fixed and of polynomial height. The reduced denominator
divides it and therefore has only \(X^{o_X(1)}\) possibilities.
The fixed-complete-value bound (GD5) gives at most
\(n^{(b-c)/4}\) outputs up to the chosen allowances, using
the lower bound \(n^c\) on the denominator radical. This
dictionary bound also holds when several tags are fixed and
the remaining slots are fixed as complete rows. In a single
query slot, the nonconstant polynomial has bounded fibres.
Since there are at most \(n^{a_i+o_X(1)}\) tags in slot \(i\),
the resulting tuple count gives the first part of (FR1).
[Fixed-complete-value bound (GD5)](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md#fixed-complete-values-and-prescribed-least-repaired-targets)

Globally, (GD2) and (GD3) on the output shell give the three
upper exponents \((b+c)/4,(b-c)/2,b/3\), up to the chosen
allowances. For the first one the pure letter has radical
cap \(c+\xi\). For the second the other two letters have
combined radical cap \(b-c\) once root sizes are binned.
[Quarter and two-letter bounds](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md#group-dictionary-count)

Apply (DR1) to the actual retained graph. Its proof uses
pairwise coprimality of the input letters, the small input tag
radical caps, and distinct input directions in each slot.
Data of degree and polynomial coefficient height fixed at
\(L\) are permitted: this changes the inner zero-set and
supported-part enumeration constants, while the number of
refinements still depends only on the input tag radical
caps and the fixed support tolerance. Thus an eligible
exponent sufficiently close to \(mr\) retains exponent
\(mr-\gamma-o_X(1)\) with at most \(n^{m\eta+o_X(1)}\)
tag vectors per output, for any fixed \(\gamma,\eta>0\).
For a fixed output and tag vector there are at most
\(n^{\sum_{j\ne i}q_j+o_X(1)}\) choices of the other complete
rows, then boundedly many choices for slot \(i\). The output
dictionary therefore has exponent at least
\(r+s-a_i-\gamma-m\eta-o_X(1)\). This gives the second
part of (FR1). The strict margin in the stated threshold
absorbs all the fixed allowances and proves the count.
[Tag recovery (DR1)](ABC_PURE_TERM_DENOMINATOR_RECOVERY.md#a-bounded-refinement-of-the-eligible-tuples)

If \(q_i\le q_*\) uniformly, (FR2) also excludes caps by a
fixed margin below
\[
 \max\{4r+2(m-2)(r-q_*)_+,\ 3r+3(m-1)(r-q_*)_+\}.
\]
On balanced input shells, (GD5) gives \(q_*\le1/6\) up to the
balance and small-budget allowances when \(L\) is asymptotic
to the input budget. For \(m=2\) and \(r\) arbitrarily close
to \(1/3\) the second threshold approaches \(9r/2\).

At \(b=4r\) in the limit of zero allowances, (FR1) forces
\(a_1=a_2=a\) for \(m=2\), and all \(a_i=0\) for \(m\ge3\).
For \(m=2\) the pure output mass is \(4a\), both other
masses are \(2(r-a)\), and \(a\le r/3\). Indeed the
pointwise and global quarter bounds must both be equalities.
Equality in the latter forces the other two masses to
coincide, since their larger one is at least half their
sum in the derivation from (GD3). The \(b/3\) exponent
gives \(a\le r/3\). For \(m\ge3\) the equations give
pure mass zero and both other masses \(2r\).

## Bivariate expansion

For \(m=2\), if \(P\) has neither of the two sum or product
composition forms in Raz–Sharir–Solymosi, Theorem 2, the
necessary cap improves to \(b\ge48r/11\). Fix both tags
in a slice containing \(n^{q_1+q_2-\kappa-o_X(1)}\)
eligible pairs by averaging over the tag vector. Write
\(Q=q_1+q_2\). The output exponent in this slice is at
least \(2Q/3\) up to allowances. If \(\max_i q_i\le2Q/3\),
use the incidence upper bound
\[
 O_{\deg P}\left((|D_1|^{2/3}|D_2|^{2/3}+|D_1|+|D_2|)|C|^{1/2}\right)
\]
on graph points, where \(D_i\) are the query sets and \(C\)
the output dictionary. If \(\max_iq_i>2Q/3\), use the bounded
fibre after fixing the smaller query coordinate. The
exceptional anchors have already been removed. Thus
\((b-c)/4\ge2(2r-s)/3\), whereas (FR1) gives
\((b+c)/4\ge r+s/2\) and \(b/3\ge r+s/2\). These imply
\(b\ge\max\{14r/3-s/3,3r+3s/2\}\ge48r/11\).
[Raz–Sharir–Solymosi, Theorem 2](https://arxiv.org/pdf/1401.7419#page=4)

If the chart has the full product output exponent \(Q\) in
such slices, the same argument gives
\(b\ge\max\{6r-s,3r+3s/2\}\ge24r/5\).
The full product count holds for multiplicatively separated
rational polynomials by clearing the fixed denominators and
using divisor-bounded factorizations of nonzero products.
It also holds for a fixed additively separated polynomial
with opposite tags under the cited affordable-modulus
hypotheses. In that case keep the fixed-coefficient and
affordability requirements of the modulus proof.
[Fixed-tag expansion](ABC_AFFORDABLE_PURE_TERM_SLICES.md#bivariate-and-trivariate-slices)
· [Opposite additive orientations](ABC_AFFORDABLE_PURE_TERM_SLICES.md#mixed-orientations-and-the-two-slot-additive-case)


## Fresh radicals in products of ratios

Return to the same input family with a fixed orientation of all three
letters \(W_i,Z_i,T_i\) in each slot, and write
\(T_i=\sigma_i W_i+\tau_i Z_i\) with fixed signs. Let
\(e_i\ne0\) be integers fixed at the chosen \(L\) and use
\[
                     \frac pq=\prod_i (Z_i/W_i)^{e_i}.
\]
The two raw pure terms \(P_1,P_2\) are products of powers of
original letters, with exactly one selected input letter
in each term per slot. Let \(A_j\) be an upper bound for
the sum in base \(n\) of the radical costs of the selected
input letters in \(P_j\), outside a prescribed initial
prime pool \(S_0\) of polynomial product height. Define
the fresh radical of \(p-q\) as the product of its primes
absent from **all three** original letters of every input row.
This includes any new primes that happen to lie in \(S_0\).

For any fixed \(\Delta>0\), tuples with \(p\ne q\) and
fresh radical at most \(n^\theta\) have a uniform
power-smaller count at sufficiently small fixed \(L\) if
\[
             0\le\theta\le 2mr-\min(A_1,A_2)-\Delta.                 \tag{FR3}
\]
The condition \(A_1+A_2\ge2mr\) is also necessary for a
near-full exponent at any bounded fresh cap \(\theta\),
up to the fixed allowances. The pairwise case
\(m=2\), \(e_1=1,e_2=-1\) is the primitive projective
secant. For balanced affordable families with arbitrarily
small record costs and \(r\) arbitrarily close to \(1/3\),
the fresh threshold approaches \(mr\), or \(2r\) for
secants. The estimates require an actual input fibre.

### Refining on all three original supports

Assume a near-full eligible exponent. Discard outputs of height
at most \(n^K\) by the tuple-chart theorem, applied to the rational
monomial, and use output height bins as above. The raw output
radical is bounded by the combined input radicals times \(n^\theta\).
This is because the two pure output letters use only input primes,
and any new prime on the third letter is accounted for by the
stipulated fresh bound. Thus at large fixed \(K\) the remaining
outputs have small relative raw budget even if the discount pool
has large radical.

For a tuple and slot \(i\), let its bad pools consist of the
other rows' prime supports and the primes of the numerator of
\[
       1-(-\sigma_i/\tau_i)^{e_i}\prod_{j\ne i}x_j^{e_j}.
\]
Discard the tuples where this expression is zero. There are
\(O_{\boldsymbol e}(N^{m-2})\) such complementary tuples per
slot, by fixing all but one of the complementary directions
and using the univariate root bound.

Start with a common pool \(S=S_0\) and fix the entire
\(S\)-supported parts of all input letters per slot at
subpower cost. For any fixed \(\eta>0\), one can refine
\(S\) and retain a nearly full exponent on which each input
letter shares radical at most \(n^\eta\) with its bad pools
outside \(S\). Indeed, if fewer than half the edges are
good, some slot-letter test has a fixed positive fraction
of bad edges. Remove vertices in that slot of degree below
a sufficiently small constant times \(E/N\). Fix a
complementary tuple with at least a constant times
\(E/N^{m-1}\) retained bad neighbours. Keep all edges incident
to those neighbours in that slot, leaving at least a constant
times \(E^2/N^m\) edges. Add the chosen bad pool to \(S\)
and fix its full supported parts in all slots. This increases
the total fixed input radical exponent by at least \(\eta\),
so there are only boundedly many steps. The enlarged pool
has polynomial product height at fixed \(L\). The fixed
parts cost an inner subpower amount at each step.
With sufficiently small initial exponent loss, the graph
therefore retains an exponent within any requested fixed
tolerance of \(mr\).
[Uniform fixed-support enumeration](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#a-radical-count-contradiction-to-the-claimed-smooth-number-lower-bound)

For a fixed reduced output, the full powers of exclusive
input free primes outside \(S\) survive in the corresponding
pure primitive terms, multiplied in valuation by \(|e_i|\).
Allocate the primes and extract their original valuations
at subpower cost: the number of slots is fixed, and the
outputs have polynomial height. More explicitly, first
record the outside-\(S\) parts of each input free letter
on the supports shared with another row. Their radicals
cost at most \(n^{2m\eta}\), and fixing their full parts
costs a subpower amount. The fixed \(S\)-parts and these
new records specify all contributions to each pure term
at the now known shared or pooled primes. The remaining
prime powers in the output must be exclusive, and each is
assigned to one of the allowed slots. Their original
valuations are then determined. The two recovered letters
give each row. The number of output values is at least
\(n^{mr-o_{\rm bud}(1)+o_X(1)}\), where the fixed positive
budget allowance \(o_{\rm bud}(1)\) can be made arbitrarily
small.

Let \(p_0\notin S\) divide the primitive difference
\(p-q\). If it belongs to an original free letter, then
it must belong to multiple row supports, since otherwise
one raw pure term is a unit and the other vanishes.
If it belongs to just one target, and no other original
letter, every input ratio is a unit at \(p_0\) and the
vanishing target gives
\(x_i\equiv-\sigma_i/\tau_i\bmod p_0\).
The vanishing difference puts \(p_0\) in the complementary
pool of that slot. This controls all original primes in the
difference outside \(S\) by the bad-pool bounds. The
outside-\(S\) radical cost of the difference is at most
\(\theta+3m\eta\), and the two pure outside-\(S\) costs
are at most \(A_1,A_2\).

### Applying the two-letter count

To use (GD3) with these partial radical costs, retain the
raw budget. Anatomic-truncate at the small relative raw
output budget. In each original root coordinate split off
the entire factor supported on \(S\), and freeze that value
at inner subpower cost using the uniform fixed-support
bound. The two resulting factors are coprime and have
the same original degree. The fixed factors have zero
charged cost and their raw caps are retained. The dyadic
caps of the remaining root coordinates sum to at most
the specified outside-\(S\) radical cost on each letter,
up to arbitrarily small truncation and shell allowances.
The ordinary original roots are squarefree and the
packed coordinate and initial coefficient errors are
absorbed by truncation. The summed original raw caps
are bounded by the small relative raw budget plus the
same allowances.

The two-letter dictionary bound now gives output exponent
at most half the minimum of
\(A_1+A_2,\ A_1+\theta,\ A_2+\theta\), up to arbitrarily
small allowances. Comparing with \(mr\) proves (FR3).
If instead an ordinary total output radical cap \(n^b\)
is given, the same tag recovery and the ordinary one-third
count force \(b\ge3mr\) for a near-full eligible monomial
graph.

These lower bounds retain the third input support through
the complementary congruence. An upper fresh radical bound
on enough actual tuples requires an additional arithmetic
input.
