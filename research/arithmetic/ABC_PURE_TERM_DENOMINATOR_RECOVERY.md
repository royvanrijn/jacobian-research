# Recovering free-letter tags from a transferred denominator

The pure-term pattern left by the
[two-mixed-factor analysis](ABC_CROSS_RATIO_FRESH_PRIME_OBSTRUCTION.md#two-mixed-factors-and-a-pure-free-letter-term)
makes the transfer a polynomial in ratios of original free
letters. A prime that occurs on just one of its input
denominators gives a pole of a controlled order unless it
divides the specialized leading coefficient. A bounded
refinement of the actual eligible tuple graph makes the
unrecovered prime support arbitrarily small. The output
denominator then recovers the input free-letter values at
low cost, even for eligible tuples of small density.

## A bounded refinement of the eligible tuples

Fix an integer \(m\ge2\), put \(n=X^L\) and work in an
affordable record of \(N=n^{r+o_X(1)}\) primitive rows
with distinct targets. Here \(r\) has a fixed positive
lower bound and is bounded above. Take the height limit
at fixed \(L\) first. Index the active slots by \(i\),
select one of the two free letters \(W_i\) in each slot,
and let \(Z_i\) denote its opposite free letter. Write
\(x_i=Z_i/W_i\). The rows give distinct directions in
each slot by primitivity. Suppose a fixed polynomial
\(P(x_1,\ldots,x_m)\) over \(\mathbb Q\) has positive
degree \(d_i\) in every slot, and that
\(\operatorname{rad}(W_i)\le n^R\) for a fixed \(R\)
on the input box. Let \(\mathcal E\) be any collection
of eligible tuples from the \(m\)-fold product.
It need not contain a product subfamily.

For every fixed \(\eta,\gamma>0\) there is
\(\kappa>0\), depending on \(m,R,\eta,\gamma\), with
the following property. Take \(\kappa\) also smaller than
half the fixed positive lower bound on \(r\). If
\(|\mathcal E|\ge n^{mr-\kappa}\), a subcollection of
size at least \(n^{mr-\gamma+o_X(1)}\) has, for each
specified rational output, at most
\[
                             n^{m\eta+o_X(1)}                 \tag{DR1}
\]
different input tag vectors \((W_1,\ldots,W_m)\).
The constants are uniform at fixed \(L\) in the
affordable modulus and the eligible set.

To prove this, denote by \(C_i(x_{-i})\) the leading
coefficient of \(P\) as a polynomial in \(x_i\), with
complementary variables denoted by \(x_{-i}\).
Discard tuples where any \(C_i=0\). The Cartesian
polynomial zero bound over the distinct rational
directions discards \(O_P(N^{m-1})\) tuples, a
power-small loss relative to \(\mathcal E\) for small
enough \(\kappa\). Put all primes of the fixed
coefficient denominators into one common pool,
and freeze full parts of the tags on that pool
at subpower cost over the edges. At each later
stage the pool remains common to all slots and
the supported parts \(H_i\) of \(W_i=H_i u_i\)
are fixed per slot. The residual \(u_i\) is prime
to the common pool.

For a fixed complementary tuple \(x_{-i}\), the bad
prime pool for slot \(i\) consists of the primes of
the other tags and of the numerator of the nonzero
\(C_i(x_{-i})\). Its prime product is of polynomial
height. A tuple is bad in slot \(i\) if the radical of
the full part of \(u_i\) supported on this pool exceeds
\(n^\eta\). If at least half of the current \(E\)
edges are good in every slot, keep them and stop.
Otherwise some slot has at least \(E/(2m)\) bad
edges. Delete vertices in that slot of degree less
than \(E/(4mN)\), losing at most \(E/(4m)\) edges.
One fixed complementary tuple now has at least
\(E/(4mN^{m-1})\) bad neighbours in slot \(i\).
Restrict to all edges incident to those neighbours.
There are at least \(E^2/(16m^2N^m)\) of them by the
retained degree cap. Record the full tag parts
supported on the selected bad pool in every slot.
This loses \(X^{o_X(1)}\) by uniform fixed-support
enumeration. Add those primes to the common pool.

Each such step increases the fixed logarithmic
radical cost in base \(n\) of the selected slot
by more than \(\eta\). The total across the slots
is at most \(mR\). Thus the number of steps is
bounded in terms of \(m,R,\eta\), independently of
\(L\). A loss from exponent \(mr\) at most doubles
at each step, up to inner subpower factors.
Taking the initial \(\kappa\) small in terms of this
bound and \(\gamma\) proves the claimed size of the
final subcollection. The pools may be chosen from
the current graph as it is restricted. No enumeration
of all possible graph refinements is needed.
[Uniform fixed-support enumeration](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#a-radical-count-contradiction-to-the-claimed-smooth-number-lower-bound)

For a prime \(p\mid u_i\) outside its bad pool,
every \(x_j\) with \(j\ne i\) is \(p\)-integral and
\(C_i(x_{-i})\) is a \(p\)-unit. All other coefficients
after the specialization are \(p\)-integral.
Primitivity gives \(v_p(x_i)=-v_p(u_i)\).
The leading term now has unique minimum valuation,
so \(v_p(P)=-d_i v_p(u_i)\). In particular the full
part of \(u_i\) outside the bad pool divides the
reduced denominator of the output. Given one output,
record the radical of the missing part of each
\(u_i\), at cost at most \(n^\eta\) per slot.
The exact powers on that support have subpower
multiplicity. The rest of each tag is a divisor
of a fixed polynomial-height integer. This proves
(DR1), with all estimates made on the actual
retained edges.

## Outputs above a fixed power of the radical scale

Return to the two-mixed-factor pattern. Call a tuple
eligible if its primitive output has radical at
most \(n^b\) with \(b\le3r\) and has nonzero terms.
The pure term is proportional to
\(\prod_i W_i^{e_i}\), with \(e_i\ge1\).
Permuting projective coordinates to make it the
denominator gives precisely the polynomial chart
in (DR1). There is a fixed \(K>0\) such that the
eligible tuples with primitive output height at
least \(n^K\) have a power-smaller count, uniformly
for sufficiently small fixed \(L\).

To see this, regularize the pure-letter degrees one
slot at a time by the original dyadic degree classes.
The product partition has a polylogarithmic number of
cells. A cell retaining an eligible exponent
\(mr-o(1)\) has a vertex class of exponent \(r-o(1)\)
in every slot, since each of the other slots has
size at most \(n^{r+o_X(1)}\). For a class with degree
\(n^{q_i+o_X(1)}\), its tags have dictionary exponent
\(a_i=r-q_i+o(1)\): the upper bound uses the number
of tags of that original degree in the whole family,
and the lower bound uses the retained class size.
Each tag has at most \(n^{q_i+o_X(1)}\) rows in that slot.
With different product cells the opposite marginal
inequality need not combine across those cells.
The proof below uses only the individual marginal
bounds and the one-slot count, so it does not require
that combination.
[Marginal regularization and the pure factor charge](ABC_AFFORDABLE_PURE_TERM_SLICES.md#dyadic-regularization-of-the-two-free-letter-marginals)

Let \(s=\sum_i a_i\). By the unary counting argument,
each tag contributes at least \(a_i-o(1)\) normalized
output radical on disjoint primes for typical tuples
of the product cell. The two mixed factors contribute
at least \(2r-o(1)\) outside a power-smaller exceptional
set of tuples on the original full Cartesian product.
Restrict this estimate to the regularized cell. A cell
retaining a nearly full eligible exponent still has
only a power-smaller exceptional subset.
Therefore a nearly full eligible exponent forces
\(b\ge2r+s-o(1)\).

Fix all pure tags and the other complete rows for a
query in slot \(i\). The raw pure term is a fixed
integer of polynomial height and its primes can be
priced at zero. On the retained typical tuples the
effective price of the output is at most \(b-s+o(1)\)
in base \(n\). The primitive pure letter has price
zero and \(b-s\ge2r-o(1)\). Binning output heights
therefore gives a dictionary cap of exponent
\((b-s)/3-\delta\) for some fixed positive \(\delta\),
using the zero-price letter gap. For a fixed output,
the query row has bounded multiplicity by the
nonvanishing of \(C_i\). Thus for each \(i\) an
eligible exponent \(mr-o(1)\) forces
\[
                       r-a_i\le (b-s)/3-\delta+o(1).
                                                                    \tag{DR2}
\]
For clarity, in these applications split
\(\log_X\) of the output height into bins of width
\(\xi L\) starting at \(KL\). There are finitely
many bins at each fixed \(L\) since the output has
polynomial height. The absolute radical cap
\(n^b\) becomes a height-relative budget at most
\(b/K\), and its effective cap is bounded below
relative to that budget in the applications of the
zero-price letter gap. Choose the fixed \(K\)
large enough for those budget thresholds, and
\(\xi\) and the inner allowances sufficiently
small. The loss \(\delta\) in base \(n\) is uniform
since \(r\) is bounded below and \(b\) bounded
above. The same binning with the ordinary
one-third bound gives an upper exponent
\(b/3\) plus any prescribed positive fixed
allowance for the global output dictionary upon
choosing \(K\) sufficiently large. Fix the
allowance smaller than the positive zero-price
gap and then fix \(K\).
[Effective one-third and zero-price letter gap](ABC_EFFECTIVE_RADICAL_BUDGET.md#a-strict-effective-saving-when-one-letter-is-free-to-record)

In particular (DR2) and \(b\le3r\) imply
\(a_i\ge\delta-o(1)\) in every slot. Apply (DR1)
to the eligible graph. For a fixed output and
input tag vector, fixing the complete rows in
all slots except \(i\) gives at most a bounded
number of completions. It follows that the
number of output values is at least
\[
        n^{mr-\gamma-m\eta-\sum_{j\ne i}(r-a_j)-o(1)}
           = n^{r+\sum_{j\ne i}a_j-\gamma-m\eta-o(1)}.
\]
Take \(\gamma,m\eta\) smaller than the uniform
positive marginal lower bound. This exceeds
\(n^{b/3+o(1)}\), a contradiction. The comparison
uses strict fixed exponent margins, so the
eligible tuple fraction is power-small. This
argument covers every fixed polynomial chart of
the pure-term pattern at these output heights.

## Small outputs for bivariate charts of the same orientation

Suppose \(P(x,y)\) has both pure tags on \(U\), with
\(x=V_1/U_1\), \(y=V_2/U_2\), and write
\(U_i=F u_i,V_i=G v_i\). Let the affordable numerical
margin of \(FG\) over the shorter free height be
at least \(\eta_0>0\). Fix any \(K>0\).
Small outputs of rational height at most \(n^K\)
have subpower fibres per fixed \(U_1\)-tag and
fixed second row, apart from the bounded set of
second-row directions making the leading coefficient
in \(x\) vanish. This holds for any such fixed
bivariate polynomial, without a special-form
hypothesis.

Clear coefficient denominators by a fixed \(J\).
If \(\log_X G\ge\eta_0/3\), both variables are
\(p\)-integral and vanish modulo \(p^{v_p(G)}\) at
every \(p\mid G\). Hence \(J P(x,y)\) has the same
value modulo \(G\) as \(J P(0,0)\). For any two
small rational outputs \(m/q\), \(m'/q'\), this
implies \(G/\gcd(G,J)\mid mq'-m'q\). The latter
integer has absolute value at most \(2n^{2K}\).
For sufficiently small \(L\) it must vanish.
There is just one small output in this case.

Otherwise the affordability margin and \(F\mid U\)
imply that \(V\) is the shorter free letter, and
\(F\) exceeds the \(v_i\)-range by a fixed numerical
margin. Fix the first tag \(U_1=A\) and the whole
second row. Write \(\Phi(Z)=P(GZ/A,y)\).
All output denominators divide a fixed
polynomial-height integer, for example
\(J A^{d_x}U_2^{d_y}\) where \(d_x,d_y\) are
the individual degrees. Record the denominator
\(q\le n^K\) at subpower cost by the divisor bound.
Relative to any one query \(v'\) in this class,
\(\Phi(v)-\Phi(v')=(m-m')/q\).

For \(p^f\Vert F\) let \(C_p\) bound the valuations
of the nonzero fixed coefficients of \(JP\), with
\(C_p=0\) outside a fixed finite set. Consider \(p\)
with \(f>C_p\) and, on that finite set, enlarge
the cutoff if necessary to absorb \(v_p(J)\).
At such a prime \(G,V_2\) are units and
\(v_p(y)=-v_p(U_2)\le-f\). In each nonzero
coefficient of \(Z^i\), the term of highest degree
in \(y\) has unique minimum valuation. Taking
\(i=d_x\) shows that the smallest valuation among
the nonconstant coefficients of \(\Phi\) is
\(-M_p\) with \(M_p\ge d_x f-O_p(1)\);
the error is zero outside a fixed finite set.
Normalize \(Q_p(Z)=p^{M_p}(\Phi(Z)-\Phi(0))\).
It has integral \(p\)-adic coefficients and a unit
nonconstant coefficient. The polynomial
\(Q_p(Z)-Q_p(v')\) therefore has Gauss norm one,
and its value at \(v\) has valuation at least
\(M_p-v_p(q)\).

An integral polynomial over \(\mathbb Q_p\) of
degree at most \(d_x\) and Gauss norm one has its
integer points of value valuation at least \(t>0\)
in at most \(d_x\) residue balls modulo
\(p^{\lceil t/d_x\rceil}\). Indeed factor it over
\(\overline{\mathbb Q}_p\). At an integral input
the nonintegral root factors contribute the same
absolute values as their roots, cancelling the
absolute value of the leading coefficient by
multiplicativity of the Gauss norm. At least one
integral root must consequently be within
\(p^{-t/d_x}\). The intersection of its ball with
\(\mathbb Z_p\), if nonempty, is one residue class
at the indicated precision. This argument also
handles repeated roots. Here it gives at most
\(d_x\) classes modulo
\(p^{(f-v_p(q)-O_p(1))_+}\). Omitted exceptional
prime powers cost only fixed factors.
The Chinese remainder theorem combines the classes
to modulus \(\gg_P F/q\), with
\(d_x^{\omega(F)}=X^{o_X(1)}\) classes. For small
fixed \(L\), this modulus still exceeds the entire
\(v\)-range. Hence there are subpower-many queries
in total for the fixed tag and anchor.

Finally, on any dyadic tag-degree class having
a near-full exponent of low-height eligible pairs,
the radical lower bound \(2r+a_1+a_2-o(1)\le b\le3r\)
forces at least one \(a_i\le r/2+o(1)\). In that
slot, the fixed-tag and fixed-anchor subpower bound
counts at most \(n^{r+a_i+o(1)}\) pairs. This is
power-smaller than \(n^{2r}\). Combining with the
large-output case excludes dense eligible pairs in
the same-orientation two-slot pure-term pattern
without a numerical lower bound on output height.
It remains a restriction on proposed transfers of
the original fibres. It supplies no upper radical
bound on the transferred outputs.


## Small outputs in additively separated charts with opposite tags

For opposite pure-letter orientations, the same small-output
conclusion holds when the two-variable polynomial has the
form \(g(\phi(x)+\psi(y))\), with nonconstant univariate
polynomials. Use the rational normalization of the
[opposite-orientation argument](ABC_AFFORDABLE_PURE_TERM_SLICES.md#mixed-orientations-and-the-two-slot-additive-case).
The inner sum \(s=\phi(x)+\psi(y)\) is rational.
For a fixed nonconstant rational polynomial \(g\),
\(H(g(s))\gg_g H(s)^{\deg g}\).
In lowest terms write \(s=a/b\). After clearing fixed
coefficient denominators, the homogeneous numerator
of \(g(a/b)\) and its denominator \(J b^{\deg g}\)
have gcd bounded in terms of \(g\): for primes of \(b\)
of sufficiently high valuation the leading term of the
numerator has its fixed leading-coefficient valuation,
and all other exceptions are at bounded powers of fixed
primes. The maximum absolute value of the two
homogeneous coordinates is \(\gg_g\max(|a|,b)^{\deg g}\)
by the dominant leading term when \(|a|/b\) is large,
and the denominator otherwise. This proves the height
inequality. Thus if the output height is at most \(n^K\),
the inner sum has height at most \(n^{K'}\) for a
fixed \(K'\).

Fix the tag in one query slot and the other complete
row. Rename the query coordinates \(x=Gv/U\), and
write the varying part of \(s\) as \(\phi(x)\).
The denominators of \(s\) divide a single polynomial-height
integer because the two tags are fixed. Record the
denominator \(q\) at subpower cost and compare to one
query \(v'\) in its class. This gives
\(\phi(Gv/U)-\phi(Gv'/U)=(m-m')/q\), with
\(|m|,|m'|,q\le n^{K'}\).
If \(G\) has positive numerical weight at least
\(\eta_0/3\), the difference of the two numerators is
divisible by \(G\) up to a fixed coefficient factor.
It vanishes for small \(L\), and \(\phi\) has bounded
fibres. If \(G\) has smaller weight, \(V\) is the
shorter free letter and \(F\) is larger than the
\(v\)-range by a fixed margin. Apply the Gauss-norm
argument of the preceding section to
\(\Phi(Z)=\phi(GZ/U)\). There is no specialized leading
coefficient from the other slot. A nonconstant
coefficient has valuation at most
\(-(\deg\phi)v_p(F)+O_p(1)\) at \(p\mid F\),
and the same modulus \(\gg F/q\) recovers \(v\)
with subpower multiplicity.

Choose a query slot with the smaller tag exponent.
The bound \(b\ge2r+a_1+a_2-o(1)\) again makes this
exponent at most \(r/2+o(1)\) for any candidate
near-full eligible exponent at \(b\le3r\).
The small-output pairs are power-small, as in the
same-orientation case. The large-output bound
above handles the other pairs. Thus every two-slot
additive composition in the pure-term pattern is
excluded from a dense eligible graph at \(b\le3r\),
without an output height hypothesis.
