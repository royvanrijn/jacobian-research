# Radical prices from binary cyclotomic factors

The factorization of a difference or sum of two exact weighted
powers exposes several separately countable integers. A nonlinear
factor can be fixed at its own radical price and has only
height-subpower many primitive root pairs. Linear factors use a
different argument. They yield additional additive comparisons
whose radical costs can be averaged. The result saves for
factor-rich binomials even without an inexpensive full divisor
on the third letter.

## Individual nonlinear factor fibres

Fix \(m\ge2\) and put \(n=X^L\).
Consider primitive positive additive triples in \(X/2<C\le X\)
with two designated letters
\[
             U=A x^m,\qquad V=B y^m .
\]
The third letter is \(W=|A x^m+\epsilon B y^m|>0\)
with a fixed sign \(\epsilon=\pm1\). Write
\(F(S,T)=A S^m+\epsilon B T^m\).
By primitivity \(A,B\) are coprime, \(x,y\) are
coprime, and \((U,V,W)\) are pairwise coprime.
At fixed \((A,B)\) write
\[
             F(S,T)=\pm\prod_{j=1}^{J}F_j(S,T)
\]
as a product of primitive irreducible integral homogeneous
forms of positive degree. All factors are distinct,
with no common root in \(\mathbb P^1(\overline{\mathbb Q})\):
the binomial has nonzero discriminant. A factor has
nonzero value at the given roots, has absolute value
at most \(X\), and divides \(W\).
The factors can be ordered by any fixed convention.
Their coefficients are at most \(O_m(X^2)\). For example,
each leading coefficient divides \(A\), and the
complex roots of \(F(S,1)\) have absolute value at
most \(X^{1/m}\). Apply the elementary-symmetric
coefficient bound to each factor.

Fix a prescribed prime pool \(\mathcal P\) of polynomial
height containing all primes dividing \(ABm\). The
evaluated factors have pairwise disjoint prime supports
outside \(\mathcal P\): a prime \(p\nmid ABm\) dividing
two factors would give a repeated projective root of
\(F\bmod p\). Indeed a primitive pair represents a point
of the projective line modulo \(p\), neither coordinate
can vanish at a root of this binomial, and the derivative
there is a unit. Their outside supports are also disjoint
from those of \(x,y\) by original primitivity. If the
coefficient pair moves, choose \(\mathcal P\) after
fixing that pair, so the supported part remains a
known record.

**BFR1.** For a nonlinear factor \(F_j\) a cap
\(\operatorname{rad}_{\mathcal P^c} F_j(x,y)\le n^{r_j}\)
gives at most \(n^{r_j}X^\varepsilon\) lifts for every
\(\varepsilon>0\), at fixed coefficients. Choose its
outside radical and enumerate the powers on that
support and on \(\mathcal P\) at height-subpower cost.
For the fixed nonzero value, a factor of degree
at least three has \(X^{o_X(1)}\) primitive roots by
Stewart's Thue bound. The form has content one and
nonzero discriminant and the bound is uniform in
its coefficients. For degree two write
\(F_j=aS^2+bST+cT^2\), with \(a\ne0\). Completing the
square gives
\[
      (2ax+by)^2+(4ac-b^2)y^2=4aF_j(x,y).
\]
The discriminant is nonzero, all coefficients and
the right side are of polynomial height, and the
map on pairs is injective. The coefficient-uniform
quadratic norm bound gives \(X^{o_X(1)}\) lifts in
that case as well. All extra count allowances are
taken at the fixed budget before the height limit.
[Stewart, Theorem 1 and its consequence on printed p. 796](https://uwaterloo.ca/pure-mathematics/sites/default/files/uploads/documents/s0894-0347-1991-1119199-x_0.pdf)
· [Browning–Verzobio, Lemma 2.3](https://arxiv.org/html/2608.24512v1#S2.SS1)
· [Prescribed-support enumeration (EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

## Averaging linear and nonlinear prices

**BFR2.** Fix a weighted common-order chart as above
with an original total raw radical cap \(n^R\), where
\(R\) is fixed, retaining any original radical or
repair constraints. Suppose the coefficient dictionary
has exponent \(f\) in base \(n\), and every retained
row obeys
\[
 \log_n\operatorname{rad}_{\mathcal P^c}(UVW)\le s,
 \qquad \log_n\operatorname{rad}_{\mathcal P^c}(W)\le t .
\]
Use lower bounds \(J_0\) for the number of irreducible
factors and \(J_{2,0}\) for the number of its nonlinear
factors on the actual coefficient states. For all
sufficiently small fixed \(L\), the count exponent
in base \(n\), up to any fixed positive allowance
and \(X^\varepsilon\), is at most
\[
        f+\min\!\left\{\frac{s}{J_0+2},
                               \frac{t}{J_{2,0}}\right\},
\]
with the second expression used when \(J_{2,0}>0\).
The fixed-order and height data and the required
count margin determine how small \(L\) must be.
A finite allowed list of common orders at each
fixed budget is permitted when the proof's stated
height and small-budget hypotheses hold for the
list used. In particular for the linear-factor derived rows,
keep \(mL\) sufficiently small compared to the fixed raw
caps and the MCB2 count tolerance. Fixed \(m\) satisfies
this when \(L\) is sufficiently small.

To prove the assertion refine all the outside
radical masses on a fixed fine mesh. The choices
are finite at the fixed order. Fix the coefficient
pair and use \(d'\) for the within-record count
exponent in a mass bin, before count allowances.
By BFR1, every nonlinear factor has outside cap
at least \(d'\). Their sum is at most \(t\),
which gives the second bound.

If there are no linear factors, the original
equation and the arbitrary-fixed-coefficient
two-cheapest bound give
\(2d'\le r_x+r_y\), with \(r_x,r_y\) the two
outside root radical caps. Fix the supported
parts on \(\mathcal P\) before applying that bound.
The complete height in the original shell is
positive relative to \(X\), and the remaining
variable raw cap is bounded by \(R\) up to
arbitrarily small allowances. Sum this inequality
with the nonlinear factor prices.

If \(F\) has a rational linear factor, then
\(-\epsilon B/A\) is an \(m\)-th power in
\(\mathbb Q\). Coprimality of \(A,B\) implies
that each is a positive integral \(m\)-th power.
Absorb their roots into \(x,y\) to write
\(U=u^m,V=v^m\). The new roots are coprime.
Their outside radicals are those of the old
roots since the coefficient primes are in the
pool. Gauss's lemma identifies the linear
factor values up to sign after this rescaling:
with one linear factor its absolute value is
one of \(|u-v|\), \(u+v\); with two factors
both occur. The two-factor case has the
difference sign and \(m\) even.

For one linear factor \(T\), the three values
\(u,v,T\) are, after permutation, an actual
positive primitive additive triple. Each of
the three choices of two charged roles gives
\(2d'\le r_i+r_j\) by the same half bound.
The largest height of the derived triple is
at least \((X/4)^{1/m}\): in the original shell
at least one of the two selected letters is
at least \(X/4\). Its value recovers both roots
up to a bounded orientation choice. Summing
the three inequalities gives
\(r_u+r_v+r_T\ge3d'\).

For two linear factors set \(T_-=|u-v|\)
and \(T_+=u+v\). Any three of
\(u,v,T_-,T_+\) give a positive addition
after moving fixed factors of two into their
roles:
\[
 T_-+\min(u,v)=\max(u,v),\quad
 u+v=T_+,\quad T_-+T_+=2\max(u,v),\quad
 T_-+2\min(u,v)=T_+ .
\]
All these rows have gcd at most two and
largest height at least a fixed multiple
of \(\max(u,v)\). They reconstruct the two
roots up to the bounded division and ordering.
Any pair of the four labels can be charged
via one of these rows. Averaging the six
half-price inequalities gives
\(r_u+r_v+r_{T_-}+r_{T_+}\ge4d'\).
Again sum the linear and nonlinear contributions.
The sum of all outside raw caps is at most
\(s\) up to the mesh loss by the disjoint
support calculation. This gives
\(d'\le s/(J+2)\) in all cases. The two caps
on derived height and on actual remaining raw
radical hold at the coefficient-dependent
support state. Refine the mesh within the
fixed allowance and add the dictionary price.
[Arbitrary-coefficient two-cheapest-role bound (MCB2)](ABC_MULTIPLE_SHORT_COFACTOR_BOUNDARY.md#bounded-group-dictionaries-with-arbitrary-fixed-coefficients)
· [Four additive identities and factor radical compression](ABC_EXACT_POWER_DIFFERENCE_RADICAL_COMPRESSION.md)

### Sharper shift-price comparisons with linear factors

In a coefficient state with exactly one linear factor
of cap \(r_T\), the two pairs \((u,T),(v,T)\) give
\[
 4d'\le r_u+r_v+2r_T\le
             s+t-2\sum_{\deg F_j\ge2}r_j .
\]
Thus a further bound is \(d'\le (s+t)/(2J+2)\).
With two linear factors, count directly with the
pair \(T_-,T_+\) using \(T_-+T_+=2\max(u,v)\).
Then \(2d'\le r_{T_-}+r_{T_+}\), and the remaining
\(J-2\) nonlinear factors each cost at least \(d'\).
It follows that \(d'\le t/J\), which also holds
when there are no linear factors.

For an odd common order divisible by three the
one-linear-factor case also gives \(d'\le t/J\).
In the difference orientation its linear factor
is \(T=|u-v|\) and \(Q=u^2+uv+v^2\) divides
the shift. The relation \(T^2+3uv=Q\) is primitive
up to a gcd dividing three. In the sum orientation
use \(T=u+v,\ Q=u^2-uv+v^2\) and \(T^2=Q+3uv\).
The roots are recovered up to bounded orientation
choices from \(T,uv\). The two roles \(T^2,Q\)
therefore give \(2d'\le r_T+r_Q\) at positive
derived height. Add the other \(J-2\) factor prices.
These refinements are assessed on the same coefficient
and mass bins as BFR2 before adding the catalogue
price \(f\).

For example, if two letters are exact \(m\)-th
powers, the cyclotomic factor count \(J\) is
\(\tau(m)\) for a difference and
\(\tau(2m)-\tau(m)\) for a sum, where \(\tau\)
counts positive divisors. Hence the coefficient
\(1/(J+2)\) holds for their entire small-budget
subcount, with arbitrary fixed count allowance
and \(\varepsilon\). This includes two square
summands, where the single quadratic factor
gives coefficient \(1/3\). A smaller shift cap or
unbalanced original caps can improve that instance. Two exact powers with
a difference at any fixed order \(m\ge2\)
give a strict improvement, with no hypothesis
on the third letter's exact valuations or on
a full divisor. In a balanced difference
chart the separate nonlinear-factor arm
has coefficient at most
\(1/[3(\tau(m)-1-\mathbf1_{2\mid m})]\)
when that denominator is positive. These
conclusions use actual factor multiplicity.
Capelli's criterion describes when reducibility
can start for a weighted binomial. Arbitrary
coefficient states can have irreducible forms,
so the larger factor count is an additional
arithmetic condition.
[Browning–Verzobio, Lemma 2.5](https://arxiv.org/html/2608.24512v1#S2.SS1)


## Fresh cyclotomic radicals on original rows

**BFR3.** Take a balanced actual primitive family in
\(X/2<C\le X\), of row exponent \(d=1/3\) in base
\(n=X^L\), with all three original caps
\(n^{d+o(1)}\). Let \(\mathcal E\) be a fixed finite
set of distinct integers at least two, and let
\(\mathcal P\) be any prescribed polynomial-height
prime pool. Enlarge \(\mathcal P\) by primes
dividing \(\operatorname{lcm}(\mathcal E)\).
For every fixed \(\gamma>0\), outside a power-small
subfamily at sufficiently small fixed budgets,
every one of the actual outputs
\[
                       \Phi_j(C,B),\qquad j\in\mathcal E ,
\]
has outside-\(\mathcal P\) radical greater than
\(n^{d-\gamma}\). These outside primes are
disjoint across \(j\) and disjoint from the support
of the three original letters. The statement is
uniform in the prescribed pool and requires no
lower bound for complete-letter multiplicity.

For \(j\ge3\) the cyclotomic form has degree
\(\varphi(j)\ge2\), content one and nonzero
discriminant. Count its fixed values at the outside
radical price and apply the nonlinear factor fibre
argument of BFR1 with \(x=C,y=B\). The ordered
pair recovers the entire original row. For
\(j=2\) use the actual positive relation
\[
                          (C-B)+(C+B)=2C .
\]
Its gcd divides two. Charge the half-price using
the original cap on \(C-B=A\) and the proposed
outside cap on \(C+B\), after fixing the
\(\mathcal P\)-supported powers. The maximum
is comparable to \(X\). This bounds the low
radical subfamily by \(n^{d-\gamma/2+o(1)}\).
For either case compare to the original row
exponent and take the chart and counting
allowances inside the power margin.

For disjointness, all these forms and \(C-B\)
are distinct factors of \(C^M-B^M\), where
\(M=\operatorname{lcm}(\mathcal E)\).
Outside primes dividing \(M\) their supports
are disjoint by squarefreeness modulo \(p\).
A prime dividing \(CB\) cannot divide any
of the evaluated factors by primitivity and
the nonzero unit end coefficients of the
cyclotomic forms. This proves the asserted
freshness, even if the original family carries
extra eligibility restrictions. In particular,
the total fresh radical across any fixed
number of distinct cyclotomic outputs exceeds
any given bounded normalized cost when that
number is sufficiently large. The number of
factors and the budget cutoffs are fixed before
the height limit. There is no bound on the
original exceptional count from this implication
alone.
[Binary cyclotomic fibres (BFR1)](#individual-nonlinear-factor-fibres)
· [Actual dyadic factor identities](ABC_EXACT_POWER_DIFFERENCE_RADICAL_COMPRESSION.md#arbitrarily-many-groups-indexed-by-powers-of-three)

## Joint fresh radicals across all three orientations

In the shell \(X/2<c\le X\) of positive primitive rows
\(a+b=c\), take \(n=X^L\), with the height limit at each
fixed \(L\). For each \(T\in\{a,b,c\}\), prescribe a fixed
finite index set \(E_T\subset\{2,3,\ldots\}\). For \(d\in E_T\)
write the corresponding output as
\[
 H_{a,d}=\Phi_d(c,b),\qquad
 H_{b,d}=\Phi_d(c,a),\qquad
 H_{c,d}=|\Phi_d(a,-b)| ,
\]
using homogeneous cyclotomic polynomials. The case \(a=b\)
gives only the primitive row \((1,1,2)\). We discard it
before defining \(H_{c,2}=|a-b|\). The common output
\[
                   H_{a,6}=H_{b,6}=H_{c,6}=a^2+ab+b^2
\]
is included only once wherever these labels occur.
Every other output retains its label.

Let \(J\) be the number of retained labels of indices at
least three and let \(I=\{T:2\in E_T\}\), with
\(\ell=|I|\). Assume \(J+\ell>0\).
Choose a prescribed prime pool \(\mathcal P\) supported on
an integer of polynomial height, uniformly in the family,
and fix the caps
\[
 \operatorname{rad}(abc)\le n^B,\qquad
 \operatorname{rad}_{\mathcal P^c}
       \left(\prod H_{T,d}\right)\le n^t
\]
with \(B,t\ge0\). Optionally bound the separate original
radical masses by \(\rho(T)=\log_n\operatorname{rad}(T)\le r_T\).

**BFR4.** At every sufficiently small fixed \(L>0\), the
count is \(O(n^{U+\eta}X^\varepsilon)\) for any prescribed
\(\eta,\varepsilon>0\), where one may take
\[
 U=\frac{t+\sum_{T\in I}r_T}{J+2\ell}.
 \tag{BFR4}
\]
With at least two linear indices, one may also use
\(U=t/(J+\ell)\). With all three, one may additionally
use \(U=(B+t)/(J+6)\).
For one orientation with any nonempty index set \(E\), one
may also use
\[
                         U=\frac{B+t}{3+|E|}.
\]
All small-budget thresholds can depend on the fixed cap
data, the indices and the allowances. The estimates are
uniform in the prescribed pool under its polynomial
height bound.

### Separation of the output supports

Put \(u=a/c\), so that \(b/c=1-u\). Up to signs, the
three outputs, divided by \(c^{\varphi(d)}\), are values
of the following integral polynomials in \(u\), respectively:
\[
 \Phi_d(1-u),\qquad \Phi_d(u),\qquad
             (1-u)^{\varphi(d)}\Phi_d\left(-\frac{u}{1-u}\right).
\]
The expressions in this display use univariate cyclotomic
polynomials. The third polynomial has nonzero value at
\(u=1\) after clearing powers of \(1-u\).
Two distinct indices in one orientation give polynomials
with no common complex root. For two different
orientations a common root would make two of
\[
                      u,qquad 1-u,qquad -u/(1-u)
\]
have complex absolute value one. This implies that all
three do and that \(u=e^{\pm i\pi/3}\). All three arguments
are then primitive sixth roots of unity. Hence the only
common factor across orientations is the shared order-six
output already identified. Every other pair of polynomials
has nonzero integer resultant.

Outside the primes dividing these fixed resultants and
the selected indices, the output supports are disjoint
from one another and from the support of the original
letters. For the latter assertion, a prime dividing a root
in a homogeneous pair cannot divide its output because the
other root is a unit and both end coefficients are units.
At a prime dividing the omitted original letter the
ratio of the two roots is one, and a prime outside \(d\)
cannot divide \(\Phi_d(1)\) for \(d\ge2\), by squarefreeness
of \(T^d-1\) modulo that prime. For the former assertion,
now with \(p\nmid c\), a common output prime gives a common
root of two of the stated polynomials at \(a/c\bmod p\).
This is excluded by the resultant. Enlarge the prescribed
pool by these fixed primes. The outside cap only decreases.
Numerical coincidences between distinct retained labels
then have their whole common support in this enlarged pool.

### Charging the selected roles

Bin the separate outside radical prices of the output
labels, and the original prices, in fixed sufficiently
fine intervals. Let \(D\) be the normalized exponent in
a bin, with the small allowances suppressed.
For each nonlinear output of outside price \(w\), the
fixed-value proof of BFR1 gives \(D\le w\). The homogeneous
form is irreducible of fixed degree and has nonzero
discriminant. Fix its outside radical and enumerate all
powers on that support and the prescribed pool. The
primitive roots recover the original row, up to the
bounded multiplicity allowed in BFR1.

For an output at index two of outside price \(w\), use
one of the genuine positive comparisons
\[
        a+(c+b)=2c,qquad b+(c+a)=2c,qquad
                   |a-b|+2\min(a,b)=c ,
\]
in the three orientations, respectively. Their common
gcd divides two, the complete height is comparable to
\(X\), and each recovers the original row up to bounded
choices, including a swap in the last comparison.
Fix all prime powers on the prescribed pool in the
three derived terms at height-subpower cost. They
enter as fixed coefficients, retaining original
primitive normalization after the bounded common
division. The half-price estimate MCB2 charges the
linear output at \(w\) and the other indicated role
at most \(\rho(T)\). The third role also has bounded
raw cap, since it is a bounded multiple of an original
letter. Thus \(2D\le w+\rho(T)\).

Sum these bounds over the selected labels. Their outside
prices sum to at most \(t\) by the disjointness just
proved. This gives (BFR4) and the \((B+t)/(J+6)\) bound
with all three linear labels. If at least two linear
outputs are selected, every pair has a comparison
with a known original third role. When \(a>b\), put
\(H_a=a+2b,\ H_b=2a+b,\ H_c=a-b\). The rows are
\[
             H_a+H_b=3c,\qquad H_c+3b=H_a,\qquad
             H_b+H_c=3a .
\]
Interchange \(a,b\) in the other ordering. Each gcd
divides three, and each pair recovers the originals
from two independent linear equations up to bounded
normalization and ordering. Each comparison has height
comparable to \(X\). The original third role has bounded
radical, so MCB2 gives \(2D\le w_i+w_j\) after
fixing all supported parts at height-subpower cost.
For two labels this is one inequality. For three,
average the three comparisons to give \(\ell D\)
at most the total linear outside price. The \(J\)
nonlinear inequalities now give \((J+\ell)D\le t\).
For a single orientation containing
index two, this sum bounds \((|E|+1)D\) by \(t+\rho(T)\).
Separately, MCB2 on the other two letters of the
original row bounds \(2D\) by their total radical cap.
Adding gives the single-orientation formula.
If the single orientation omits index two, the nonlinear
sum bounds \(|E|D\) by \(t\), and the small-budget
one-third estimate on the original count bounds
\(3D\) by \(B\) up to its arbitrarily small
normalized allowance. In the shell use relative
budget slightly above \(BL\). Taking the finite
bin widths and count allowances inside \(\eta\)
gives the stated estimates without exchanging the
height and budget limits.

For example, with all three index sets consisting
of the divisors \(d\ge2\) of a fixed \(M\) divisible
by two, denote \(K=\tau(M)\). The denominator
\(J+6\) is \(3K-2\mathbf 1_{6\mid M}\). On a putative
unrestricted saturation sequence with \(B=1\)
and normalized exponent tending to \(1/3\), any
fixed cap below
\[
                         K-1-\frac23\mathbf 1_{6\mid M}
\]
on the outside radical of this union gives a strict
saving on the capped family. These output radicals
are on fresh primes after the bounded exclusions
above. A bound on their size does not follow from
the original raw radical budget.
[The half-price estimate (MCB2)](ABC_MULTIPLE_SHORT_COFACTOR_BOUNDARY.md#bounded-group-dictionaries-with-arbitrary-fixed-coefficients)
· [Small-budget one-third estimate](ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13)
· [The original one-third limiting count](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#the-unrestricted-one-third-endpoint)


## The remaining fresh-radical estimate

For a concrete choice, take \(E_a=E_b=E_c=\{2,3,6\}\).
There are four distinct nonlinear factors and three linear
factors. The coefficient in BFR4 is at most
\(\min\{t/7,(1+t)/10\}\) at original total budget one.
Consequently a fixed loss \(\sigma>0\) below the
fresh-union radical exponent \(7/3\) gives an
original exponent of at most \(1/3-\sigma/7\),
up to any fixed positive allowance.

The missing selection statement can use any single
fixed index catalogue with \(\ell\ge2\), and set
\(q=J+\ell\). Suppose there is a sequence of fixed
budgets \(L_i\downarrow0\) along which original rows
in height shells reach exponent \(1/3-o_i(1)\) in
base \(n=X^{L_i}\) at arbitrarily large inner
heights. A sufficient selection statement asserts
that some such fixed catalogue and some fixed
\(\sigma>0\) retain exponent \(1/3-o_i(1)\) on a
subsequence of those budgets and heights while the
outside union radical is at most \(n^{q/3-\sigma}\).
The pool may be prescribed at each inner height
under a common polynomial height bound. If using a
number of prescribed pools of exponent \(b\), take
\(b<\sigma/q\) by a fixed margin. The original
three radical caps are retained. With these
quantifiers BFR4 contradicts the selection. The
statement gives a route for every \(\alpha\in[1/6,1/3]\),
including occurrence, because the count has no
positive complete-fibre or dictionary hypothesis.
The existing prime-support and residue estimates
only give lower prices for the cyclotomic outputs.
They do not supply this upper price for actual rows.

## Exact-power sectors with arbitrary common order

Fix \(B=1\) for the original total radical budget and
\(n=X^L\) in the shell \(X/2<c\le X\).
Consider primitive triples with two original letters
that are exact \(m\)-th powers for some integer
\(m\ge2\). The order and the two chosen letters can
vary by row. There are absolute \(\delta,L_0>0\) such
that for every fixed \(0<L<L_0\) this sector has count
\[
                    O_{L,\varepsilon}
                       \bigl(n^{\,1/3-\delta}X^\varepsilon\bigr).
\]
This also covers two exact powers whose literal
orders have a common divisor greater than one.
For two summands with a common order that is a power
of two, both summands are squares. The same-side
\(q=2,f=0\) case of PR7 gives a fixed strict saving.
The following argument covers the target–summand
orientation and all remaining same-side orders.
Taking the smaller margin gives the asserted bound.
[The same-side quadratic case (PR7)](ABC_PROJECTIVE_RAY_CLASS_NORM_MOMENTS.md#weighted-powers-and-the-balanced-endpoint)

Here is the division of order ranges. Fix a small
count allowance \(\eta>0\), and choose \(\kappa>0\) small
in terms of the MCB2 thresholds for this allowance and
bounded raw caps. At \(mL\le\kappa\), apply the proof of
BFR2 with coefficients one and the actual order \(m\).
For a difference the number of rational factors is
\(J=\tau(m)\ge2\). For a sum it is
\(J=\tau(2m)-\tau(m)\ge2\) when \(m\) has an odd
divisor greater than one. The only steps in BFR2
requiring small budget are the MCB2 counts. On the
original row this uses \(L\). If there are linear
factors, apply MCB2 to the derived additive comparisons
in height scale \(Y=X^{1/m}\) and budget \(mL\).
All these comparisons have height bounded above and
below by fixed multiples of \(Y\) (their largest
original powered operand is at least \(X/4\)), and
their raw caps in the unchanged base \(n=Y^{mL}\)
are bounded by a fixed constant. Hence the MCB2
threshold on this scale is uniform in \(m\).
The fixed-value factor counts and finite mass bins
are used only at each order \(m\le\kappa/L\) before
the height limit. BFR2 then gives \(n^{1/4+\eta}X^\varepsilon\)
after summing these finitely many orders at fixed \(L\).

Fix \(K>8\). For \(\kappa\le mL\le K\), the
two large-step valuation spectra sector applies to
the two powered letters. The endpoint contradiction
is uniform in the choices of the orders. If there
were no uniform saving, take a sequence of fixed
budgets tending to zero. At each budget the finite
list \(\kappa/L\le m\le K/L\) gives one order with
a near-saturating inner exponent, and the two
spectra then have uniformly positive scaled mesh.
The large-step result excludes that sequence.
Finally, for \(mL>K\), each of the two roots is at
most \(n^{1/K}\). Some root is at least two, apart
from the bounded row \((1,1,2)\), and hence
\(m\le\log_2 X\). Root enumeration gives
\(n^{2/K}X^\varepsilon\) after summing these
orders. A dyadic shell summation preserves the
strict margin in the bound.
[The two large-step spectra sector](POWERFUL_TRIPLES_ABC_ENDPOINT_TRANSFER.md#large-step-valuation-spectra-and-the-small-budget-endpoint)

## Uniform saving for rationally reducible weighted binomials

Suppose instead that two letters of a primitive row
are \(U=A x^m,V=B y^m\) for some \(m\ge2\), with any
positive integral coefficients \(A,B\), and that
\(A S^m+\epsilon B T^m\) is reducible over \(\mathbb Q\)
in the actual sign orientation. Allow the order and
both coefficients to depend on the row. The same
conclusion of a fixed saving below \(1/3\) holds,
with no bound on the coefficient dictionary.

Primitivity implies \((A,B)=1\). Dehomogenizing,
Capelli's necessary condition gives either that
\(-\epsilon B/A\) is a rational \(p\)-th power for a
prime \(p\mid m\), or that \(4\mid m\) and
\(-\epsilon B/A=-4r^4\) for a nonzero rational \(r\).
In the first case \(B/A\) is a positive rational
\(p\)-th power. Coprimality makes both \(A\) and
\(B\) positive integral \(p\)-th powers. Thus the
two actual letters are exact \(p\)-th powers.
In the sum orientation \(p\) is odd because
\(-B/A<0\). The adaptive common-order result
therefore applies in either sign orientation.

The second case is a sum. At odd primes the
coefficient exponents are divisible by four.
At \(2\), the valuation of \(B/A\) is \(2\)
modulo four. Coprimality then puts one coefficient
in the form \(4u^4\) and the other in the form
\(v^4\). Absorb the actual integral roots and
\(m/4\) into \(u,v\); the two original letters
are \(4y_0^4,x_0^4\), in some order. The fixed
binomial factors as
\[
   x_0^4+4y_0^4
      =(x_0^2-2x_0y_0+2y_0^2)(x_0^2+2x_0y_0+2y_0^2).
\]
The quadratics are distinct and irreducible,
so BFR2 with fixed coefficients gives exponent
\(1/4+\eta\) for arbitrarily small fixed \(\eta\)
at sufficiently small budgets. This exhausts
the rationally reducible sector.
[Capelli's criterion, Browning–Verzobio Lemma 2.5](https://arxiv.org/html/2608.24512v1#S2.SS1)

## A two-cap obstruction to any fixed fresh-output catalogue

The following construction uses two original \(n^{1/3}\)
caps and gives actual primitive additive rows, with
both capped valuation laws wholly diffuse. It applies
to every finite fixed collection of the three
orientations above and even allows a prescribed
catalogue of pools. In this construction the original
third radical exceeds every fixed normalized exponent
as the budget decreases.

Fix finite index sets for the three orientations, a
number \(0\le\varrho<1/3\), and \(\eta>0\) with
\(\varrho+\eta<1/3\). At each sufficiently small fixed
\(L\) and height \(X\), let a prescribed catalogue
contain at most \(n^{\varrho+o_X(1)}\) prime pools, with
a uniform polynomial bound on their support products.
Here \(n=X^L\), and support products are taken only
over distinct primes. There is a family of
\(n^{1/3+o_X(1)}\) actual primitive rows \(a+b=c\),
forming a matching in all three complete letters, with
\(\operatorname{rad}(a),\operatorname{rad}(c)\le n^{1/3}\).
For every row and every pool in the catalogue, each
prescribed nonlinear output has outside radical at
least \(n^{2/3-\varrho-\eta}\), and each prescribed
linear output at least \(n^{1/3-\varrho-\eta}\).
Thus for \(J\) distinct nonlinear labels and \(\ell\)
linear ones, the joint outside radical has exponent
at least
\[
                  \frac{2J+\ell}{3}
                      -(J+\ell)(\varrho+\eta)-o_X(1).
\]
The shared order-six label is counted once. For fixed
\(J+\ell>0\) and any fixed saving \(\sigma>0\) with
\((J+\ell)\varrho<J/3+\sigma\), choose \(\eta\) small
enough to make the required fresh ceiling
\((J+\ell)/3-\sigma\) fail on every row and for every
allowed pool.

To construct the rows, use the two full prime-palette
dictionaries in the diffuse two-cap construction, with
the same mesh and the following choice of common order.
Take a prime \(q=q_L\to\infty\) bounded above by
\(\log(1/L)\), set \(s=\lfloor\log(1/L)/(2\log q)\rfloor\),
and use \(e=q^s\) instead of the dyadic order. Then
\(s\to\infty\), \(eL\to0\), and \(e\le L^{-1/2}\).
More explicitly put \(m=\lfloor L^{-1/4}\rfloor\),
\(R=n^{1/(3m)}\),
\[
 k_i=\left\lfloor\frac{2+(2i-1)/m}{eL}\right\rfloor,
 \qquad
 H_L=\frac1{3m}\sum_{i=1}^m eL k_i .
\]
Give each letter \(m\) disjoint palettes of
\(P\asymp_L R/(\log R)^2\) primes, in \([R/8,R/4]\)
for \(a\) and in \([R/2,R]\) for \(c\). Each word
is an actual product of one prime from each palette
to exponent \(e k_i\). The dictionary sizes are both
\(M=P^m=n^{1/3+o_X(1)}\). The palette separation gives
\(c>2a\) for every pair, and all rows
\((a,c-a,c)\) are primitive with all heights
\(\asymp_L X^{H_L}\). Also \(H_L\to1\). On each
capped letter the scaled-degree measure tends
in the two nested limits to the uniform measure
of total mass \(1/3\) on \([2,4]\).
One may select a most populated dyadic height shell
at the end, changing the count by a factor bounded
at fixed \(L\).
[Diffuse power dictionaries](ABC_DIFFUSE_TWO_CAP_PRIMITIVE_SECANT_OBSTRUCTION.md#diffuse-power-dictionaries)

For one fixed pool, counting failure of the output lower
bound is checked on the full Cartesian product of the
two dictionaries. The number of polynomial-height
values with outside radical at most \(n^t\) is
\(n^tX^{o_X(1)}\), uniformly in the pool, by
prescribed-support enumeration. A nonlinear output
has only \(X^{o_X(1)}\) pairs per exact value by the
binary fixed-value bound. Its two arguments are
a unimodular linear transform of \((a,c)\), and in
particular are coprime. A linear output allows
at most \(O(M)\) pairs per value by first fixing \(a\)
and then solving the linear equation for \(c\).
Accordingly, over all pools and prescribed labels,
the number of Cartesian pairs failing any of the
desired individual lower bounds is
\(n^{2/3-\eta+o_X(1)}=o(M^2)\).
[BFR1](#individual-nonlinear-factor-fibres)
· [Prescribed-support enumeration](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

One can demand a growing normalized radical exponent
on the original third letter and on each selected
output \(\Phi_d(c,a)\) in the orientation omitting
that letter. Write \(a=u^e,c=v^e\).
The third letter \(b=v^e-u^e\) contains the \(s\)
nonlinear factors \(\Phi_{q^j}(v,u)\) for \(1\le j\le s\)
at sufficiently small \(L\), when \(q\) is odd.
For each fixed \(d\ge2\), eventually \(q\nmid d\), and
\[
       \Phi_d(v^e,u^e)
                  =\prod_{j=0}^{s}\Phi_{d q^j}(v,u).
\]
This follows by taking orders of the complex roots
under the \(e\)-th power map and then homogenizing.
At least \(s\) factors are nonlinear. At each fixed
\(L\) the factors used in the two factorizations are
a finite list. Their distinct prime supports
within each factorization are disjoint away from
the finitely many primes of the relevant index
and its fixed resultants. Impose the same individual
lower bound \(n^{2/3-\varrho-\eta}\) on every such
nonlinear factor, using the binary fixed-value
count for the primitive root pair \((v,u)\).
This discards only \(n^{2/3-\eta+o_X(1)}\) pairs
at each fixed \(L\). It gives the stated growing
outside radical, of exponent at least
\(s(2/3-\varrho-\eta)-o_X(1)\) in those orientations
and on \(b\). All finite-prime exclusions have
zero normalized mass in the inner limit.

Finally, thin the surviving pairs independently
with probability \(1/(M\log X)\). The expected
size is \((1-o(1))M/\log X\), and its variance
is at most its mean. The expected number of
collisions in either \(a\) or \(c\) is
\(O(M/\log^2 X)\). For any fixed \(b\), a single
nonlinear factor \(\Phi_q(v,u)\) must be one of
its divisors. The divisor bound and binary
fixed-value estimate therefore give at most
\(X^{o_X(1)}\) Cartesian pairs for that \(b\).
The expected number of collisions in \(b\)
is \(X^{o_X(1)}\), also negligible compared with
the expected size. Delete one pair per collision.
The resulting matching has \(n^{1/3+o_X(1)}\) rows
and retains the outside radical inequalities.
The \(o_X(1)\) errors are always taken with
\(L,q,s,m\), the cap margins, and the selected
finite indices fixed. This example gives an
obstruction with only the two original caps,
and the actual third-letter lower bound places
it outside the original total small-radical budget.
