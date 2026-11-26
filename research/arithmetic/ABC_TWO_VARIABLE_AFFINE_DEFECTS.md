# Affine defects on two-variable additive families

Fixing a complete letter makes every affine form in the other
two coordinates a constant or a one-variable affine form. A
different fixed-letter slice can test each pair of forms.
The three slice counts together give simultaneous defects
outside one query-dependent pool, including for forms in
different additive directions.

## Role-capacity bounds for affine forms

Fix \(B\ge1\), \(d>0\), and \(0<\eta<1\). Use \(n=X^L\)
with \(0<L\le1\) fixed in the height limit. Let
\(\mathcal E\) be a family of at least \(n^d\) actual positive
primitive triples \(a+b=c\le X\), with complete dictionary
sizes at most \(n^{h_T}\) for the roles \(T=a,b,c\). Write
\[
     F_{\max}=\max_{T,t}\#\{(a,b,c)\in\mathcal E:T=t\}.
\]
Additional restrictions on the original triples and their
three radicals may be retained.

Set \(B_*=B+3\) and \(\tau=\eta/(128B_*)\). Let
\(2\le m\le(\log X)^{\tau L}\) and let \(F_1,\ldots,F_m\)
be pairwise nonproportional nonconstant integer affine forms
in \(a,b\), with coefficients of magnitude at most \(X^B\).
Suppose all their values are nonzero on \(\mathcal E\).
Prescribe a prime pool \(\mathcal P_0\) of cardinality at
most \((\log X)^{1+2\tau L}\). At a query \(z=(a,b,c)\),
define
\[
 \mathcal P(z)=\mathcal P_0\ \cup\
   \{p:p\mid F_i(a,b),\ p\mid F_j(a,b)\text{ for some }i\ne j\}.
                                                               \tag{DA1}
\]
The primes in this pool can vary with the query.

Say a form *moves on a \(T\)-slice* if it is nonconstant on
the affine line \(T=t\). This condition is independent of \(t\).
Every nonconstant form moves on at least two of the three
types of slice, so each pair moves together on at least one.
Put
\[
 \begin{split}
 u_i&=\max\bigl(0,\ \max_{T:\,i\ {\rm moves}}
                                   (d-h_T-2\eta)\bigr),\\
 u_{ij}&=\max\bigl(0,\ \max_{T:\,i,j\ {\rm move}}
                                   (d-h_T-2\eta)\bigr).
 \end{split}
\]
Except for
\[
             O_{B,d,L,\eta,\{h_T\}}
                   \bigl(n^{d-3\eta/2}+m^2 F_{\max}\bigr)       \tag{DA2}
\]
queries, all the following bounds hold simultaneously:
\[
 \begin{split}
 \log_n\mathcal R_{\mathcal P(z)}(F_i)&\ge u_i,\\
 \log_n\mathcal D_{q,\mathcal P(z)}(F_i/F_j)&\ge u_{ij}
                                       &&(q\ge2,\ i\ne j),\\
 \log_n\mathcal D_{k,\mathcal P(z)}(F_i)&\ge u_i
                                       &&(k\ge 8B_*/(\eta L)).
 \end{split}                                                   \tag{DA3}
\]
Here \(\mathcal R,\mathcal D\) are the outside radical and
outside power-defect products in (AF1). The inequalities are
strict whenever the corresponding floor is positive. For
example, the second bound uses the smallest complete
dictionary among the roles on which both forms move.
[Growing one-variable affine defect count (AF1)](ABC_GROWING_AFFINE_DEFECT_FIBRES.md#growing-families-of-affine-forms)

### Slice pools and degeneracies

On an \(a\)-slice use \(b\) as the free variable, on a \(b\)-slice
use \(a\), and on a \(c\)-slice substitute \(b=c-a\) and use
\(a\). For any pair of moving forms, its determinant as a
function of the fixed value \(t\) is a nonzero polynomial
of degree at most one. Indeed the slopes in the free variable
are fixed and nonzero. An identically zero determinant would
give a fixed proportionality of their coefficients in both
the free variable and \(t\), hence a proportionality of the
original two-variable forms. Therefore at most \(O(m^2)\)
fixed values, across the three roles, give any proportional
pair of moving forms. Discarding these slices costs
\(O(m^2 F_{\max})\) queries.

On each remaining slice put the prime support of every constant
form on that slice into the prescribed pool for the one-variable
theorem. These forms have fixed nonzero values whenever the
slice meets \(\mathcal E\). Their total prime support has size
\(O_B(m\log X/\log\log X)\), since all form values have magnitude
at most \(X^{B_*}\). Adding these primes and \(\mathcal P_0\)
fits the prescribed-pool cardinality in (AF1) with parameter
\(\eta/(64B_*)=2\tau\) for sufficiently large \(X\). All slice
coefficients also have magnitude at most \(X^{B_*}\).
Adjoin the prime supports of the fixed determinants between
the moving forms, as required by (AF1).

The resulting slice pool contains \(\mathcal P(z)\) for
every query on the slice. A prime shared by two moving forms
divides their determinant. If either form is constant, the
prime lies in the full support of a constant form. Thus passing
from a slice pool to \(\mathcal P(z)\) can only increase the
radical and defect products.

In a role with \(d-h_T>2\eta\), apply (AF1) on every such slice
at \(r_T=d-h_T-\eta\), with its allowance \(\eta\).
For a slice with only one moving form, the unary proofs of
(AF1) apply directly. A slice with none supplies no test.
Each count uses the free coordinate as an injective
parameter on the actual triples of its slice. The exceptional
count summed over the at most \(n^{h_T}\) slices is
\(O(n^{h_T+r_T-\eta/2})=O(n^{d-3\eta/2})\).
Take the union over the three roles. For each positive floor
in (DA3) choose an eligible orientation attaining it. This
proves the simultaneous assertion. In particular the bound
allows the pair, order and weights to be chosen after the
query, although the list of forms is fixed before counting.

## Weighted products and radical extensions

For an integer \(q\ge2\) and arbitrary integer weights \(w_i\),
put \(k_i=q/\gcd(q,w_i)\), with \(\gcd(q,0)=q\).
Use an edge \(i,j\) of weight \(u_{ij}\) if
\(\gcd(k_i,k_j)>1\), and allow a singleton of weight \(u_i\)
if \(k_i\ge8B_*/(\eta L)\). Let \(M^*\) be the maximum
weighted sum of nonnegative edge and singleton multiplicities
with incident sum at most one at each index. On the retained
queries,
\[
 \log_n\mathcal D_{q,\mathcal P(z)}
            \left(\prod_iF_i^{w_i}\right)\ge M^*,\qquad
 \log_n\mathcal R_{\mathcal P(z)}
            \left(\prod_iF_i\right)\ge\sum_i u_i.              \tag{DA4}
\]
Indeed off \(\mathcal P(z)\) the different forms have disjoint
supports. The defect contribution from form \(i\) is
\(\mathcal D_{k_i,\mathcal P(z)}(F_i)\), with the contribution
zero in logarithmic units when \(k_i=1\). On a pair edge, choose
a prime dividing \(\gcd(k_i,k_j)\) and use the corresponding
ratio defect in (DA3). On a singleton use the unary high-order
bound. Each individual contribution is charged at most once.

If \(u=\min_{i\ne j}u_{ij}>0\), the same retained queries
satisfy the Kummer-degree lower bound
\[
 \left[K(\mu_q,F_1^{1/q},\ldots,F_m^{1/q}):K(\mu_q)\right]
                 \ge q^{\max(0,m-C)},\qquad
 C=\max(1,\lceil 2H/(Lu)\rceil),                               \tag{DA5}
\]
whenever \(K\) is a number field and a positive integer
\(J\le X^H\) is divisible by \(q|\operatorname{Disc}(K)|\).
These data may depend on the query. For each \(\ell\mid q\) prime,
two indices without a private \(\ell\)-defect prime outside \(J\)
use at least \(u\) combined normalized radical mass on \(J\).
Summing the pair inequalities over the lost indices, whose
supports off \(\mathcal P(z)\) are disjoint, gives at most \(C\)
lost indices. At retained rational primes both \(K\) and the
cyclotomic field are unramified, so normalized local valuations
in their compositum equal the rational ones. For each
\(\ell^b\parallel q\), at least \(m-C\) disjoint valuation vectors
therefore have \(\ell\)-primary order \(\ell^b\).
Kummer theory gives (DA5) via the valuation image of the
generated power classes.

With \(q\mid R\), the same bound holds for the least
\(R\)-full multipliers of \(|F_i|\) when
\(u>2B_* /(RL)\), replacing \(u\) in \(C\) by
\(u-2B_* /(RL)\). Valuations at least \(R\) account for at
most \(B_* /(RL)\) normalized radical mass per form. At each
other retained prime dividing a form the multiplier valuation
is \(R-v_p(F_i)\), preserving its nonzero defects modulo
each \(\ell\mid q\).
[Weighted defects and Kummer-degree argument (AF3, AF5)](ABC_GROWING_AFFINE_DEFECT_FIBRES.md#degree-of-the-simultaneous-radical-extension)

## Consequences for rooted completions

At the minimum complete-label scale, take
\(d=1/3-\kappa\) and \(h_a=h_b=h_c=1/6+\kappa\).
All floors in (DA3) are \(1/6-2\kappa-2\eta\) when
this number is positive, and \(F_{\max}\le n^{1/6+\kappa}\).
Thus for every fixed positive tolerance, sufficiently
small \(\kappa,\eta\) give floors above \(1/6\) minus
that tolerance on all but a power-smaller part of the
family. Zero evaluations may also be removed here.
For each form, slice by a role on which it moves and
count at most one zero per fixed value, for a total
of at most \(m n^{1/6+\kappa}\) rows.

For a fixed anchor \((a_0,b_0,c_0)\), the rooted rectangle
count (BR1) supplies at least \(n^{1/3-O(\kappa)}\)
completions \((a,b,c)\) with all four required cells.
Apply the theorem at that exponent with the same three
dictionary caps. Include the three original forms \(a,b,c\)
in the list of tests. The outside primes on each other
form are then fresh relative to all three complete query
values as well as all other forms and the prescribed
anchor pool. For the three shifted forms \(a-a_0,b-b_0,c-c_0\),
the primitive signed difference triple obtained by
dividing their common gcd retains these private packets.
Its three radicals therefore each have fresh mass at
least \(1/6-O(\kappa)-2\eta\). The ratio defects of the
primitive differences are the same as before the division.

One may also use the actual shifted outputs \(a+b_0\)
and \(a_0+b\), or \(a+b_j\) and \(a_i+b\) whenever the
specified cells are present. The anchors can be actual
capped graph labels, and every shifted output that is an
original cell retains its original radical cap. The theorem
does not infer adjacency to a growing fixed list of anchors.
Nor does the private-radical lower bound give the smaller
upper cap on the shifted product needed to contradict it.
[Rooted rectangle completions](ABC_ROOTED_PRIVATE_RECTANGLES.md#how-many-rectangles-contain-an-anchor)

## A stronger radical floor for primitive directions

In the near-minimum complete-label setting, fix a role \(U\)
and suppose an opposite letter has positive nonatomic mass
in a compact scaled-degree band with lower endpoint above
some fixed \(s\in(3,6)\). Equivalently, choose
\(3<s<\min(6,s_*)\) in the notation of (AC1).
Work on the regularized graph for that estimate and put
\(\gamma=1/6\). For a fixed actual anchor triple
\(E_0=(a_0,b_0,c_0)\), let \(E-E_0\) be the signed
coordinate difference and divide its three coordinates by
their positive common gcd \(g\). A difference with no zero
coordinate defines a primitive positive additive triple
after taking absolute values and placing the largest term
last. Denote that triple by \(D(E,E_0)\).

Fix \(b_0'>0,B'>0\) and a prescribed integer \(P\le X^{B'}\).
The anchor and \(P\) may vary with height, with \(P\) fixed
within the query family. Suppose a query class of full
limiting exponent \(2\gamma\) has
\[
    \operatorname{rad}(D(E,E_0))\le n^{b_0'},\qquad
    \mathcal R_{\{p:p\mid P\}}(D_1D_2D_3)\le n^b
\]
for fixed \(b\ge0\). Its required outside radical budget is
\[
                            b\ \ge\ 3/2-3/s .              \tag{DA6}
\]
In particular the lower bound exceeds \(1/2\) by a fixed
amount whenever \(s>3\). No bound on the raw radicals of
unrestricted directions is assumed.

Here is the proof. Fix small positive numerical margins
\(\rho,\theta\) with \(\rho+\theta<\gamma s/2\). Take all
edge, dictionary and progression-count tolerances sufficiently
small compared with these margins. The minimum-mass height
conclusion and the uniform single-value strip estimate used
for (FG3) discard \(|U-U_0|\le X^{1-\theta}\) at a fixed
power loss in base \(n\). Indeed \(U_0\ge X^{1-\theta/2}\)
for sufficiently late budgets, and the strip contains at
most \(n^{\gamma-\delta_\theta}\) distinct complete values
for some \(\delta_\theta>0\), up to arbitrarily small allowances.
Sum over a second complete dictionary to bound its edge count.
Discard zero differences in the other coordinates at the
trivial single-letter fibre bound \(n^{\gamma+o(1)}\).

For either other role \(V\), fix \(V\ne V_0\). The common
divisor \(g\) is a divisor of \(|V-V_0|\), giving only
height-subpower many choices. For each \(g\ge
X^{1-\gamma s+\rho}\), the \(U\)-values lie in
\(U_0+g\mathbb Z\) within \([1,X]\).
Equation (AC1) counts at most \(n^{\gamma-\rho/s}\)
such values up to arbitrarily small fixed exponent losses.
Summing over \(V\) and the divisors discards this branch at
a fixed power loss. The remaining direction therefore has
\[
             |U-U_0|/g\ge X^{\gamma s-\rho-\theta}.
\]
This gives a positive numerical height lower bound for the
primitive triple, independent of \(L\) along the chosen
small-budget sequence.
[Single-value strips (FG3)](ABC_DIFFUSE_FACTOR_GRAPH.md#counting-distinct-terms-and-ratios-in-strips) · [All-scale progression count (AC1)](ABC_DIFFUSE_FACTOR_GRAPH.md#an-all-scale-progression-bound-for-complete-letters)

Fix one primitive signed direction. Its original \(U\)-values
lie in a progression of step \(|U-U_0|/g\), which now has
parameter length at most \(X^{1-\gamma s+\rho+\theta}\).
Each \(U\) determines \(g\) and the whole query on that line.
A given positive primitive additive triple is the normalized
image of at most twelve signed coordinate directions.
Thus (AC1) bounds the preimage count per positive triple
by \(n^{1/s-\gamma+(\rho+\theta)/s}\) up to arbitrarily
small losses. The surviving full-exponent class therefore
produces primitive triples of count exponent at least
\(3\gamma-1/s-(\rho+\theta)/s\).

To count the derived triples from above, split their heights
into shells \(Y/2<\max_i D_i\le Y\). Here the limiting values
of \(\log_XY\) lie in a fixed compact subset of \((0,1]\)
by the numerical lower bound. Price the primes of \(P\)
at zero and all others at \(\log_Y p\). The pool is of
polynomial height in \(Y\) uniformly since \(Y\ge
X^{\gamma s-\rho-\theta}\). For each fixed \(L\), round the
relative raw budgets \(b_0'L\log X/\log(Y/2)\), and the
effective caps \(bL\log X/\log Y\), up on meshes of arbitrarily
small fixed multiples of \(L\). These meshes have finitely
many possible budgets at that fixed \(L\). Equation (EB4)
now gives at most \(n^{b/3}\) derived triples up to arbitrarily
small fixed allowances, since the height factor cancels
the normalization of the effective budget. Here the raw
relative budgets are \(O(L)\), and the shell and height
allowances are taken before \(L\to0\). Comparing bounds and
taking \(\rho,\theta\) arbitrarily small gives (DA6).
This charges primes outside a fixed anchor pool. It does not
discount supports that vary with the individual query.
[Effective radical count and fixed pools (EB4, EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)
