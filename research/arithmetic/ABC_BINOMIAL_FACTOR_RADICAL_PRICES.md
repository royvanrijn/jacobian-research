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
