# Strong Factorial from separating shared rows

Internal research proof, 21 November 2026. Repository base:
`19eb9912fb60142f8b620dac24896931d4ff3e46`.

## Result and scope

The private variables in the [balanced amplification theorem](SFC_BALANCED_PRIVATE_AMPLIFICATION.md)
can be replaced by shared groups of variables. An exact fractional packing
of separating groups supplies the required amplification. This proves a
class of actual integer row identifications, at every intermediate stage,
for arbitrary support size, arbitrary coefficients, and every moment order.

A second theorem applies directly to balanced exponent matrices whose
zero/nonzero columns are distinct. If their positive entries are large
enough, no additional variables are required. Neither theorem proves
unrestricted row identification or removal of all decorations. In
particular, equal endpoint profiles remain a hypothesis. Full SFC is open.
No novelty or external-review claim is made.

For a characteristic-zero field `K`, let `Lambda(X^a)=prod_i a_i!`.
A nonzero polynomial with `r` monomials satisfies SFC if every window
`Lambda(F^n),...,Lambda(F^(n+r-1))`, `n>=1`, contains a nonzero entry.
All inequalities below concern positive rational structure constants;
they do not impose positivity on `K` or on the coefficients of `F`.

## 1. A separating-layer certificate

Let `M_1,...,M_s` be monomials in base variables `Y`. Assume their exponent
columns have the same multiset of positive entries. The monomials may
coincide. For each layer `l=1,...,L`, choose a map

\[
 \sigma_l:\{1,\ldots,s\}\longrightarrow A_l
\]

to a finite alphabet and an integer `D_l>=1`. Introduce distinct variables
`Z_(l,a)`, independent of the base variables, and put

\[
 F(c)=\sum_{j=1}^s c_j M_j(Y)
                 \prod_{l=1}^L Z_{l,\sigma_l(j)}^{D_l}.                 \tag{1}
\]

A set `B` of layers is **separating** when the map
`j -> (sigma_l(j))_(l in B)` is injective. A certificate consists of
finitely many separating sets `B` and nonnegative rational numbers `w_B`
such that

\[
 \sum_{B\ni l}w_B\le D_l\quad(1\le l\le L),\qquad
 W:=\sum_B w_B\ge6s.                                                   \tag{2}
\]

**Theorem 1.** Under these assumptions every nonzero specialization of
(1) satisfies SFC over every characteristic-zero field, with the actual
number of nonzero coefficients as its window length.

The certificate is a finite set of rational inequalities and injectivity
checks. The theorem quantifies over all moment orders; its proof is not a
bounded search or a numerical continuation argument.

### 1.1. The counting inequality

For a composition `k=(k_1,...,k_s)` of `m`, define

\[
 C_k=\frac{m!}{\prod_j k_j!},\qquad
 K_{l,a}=\sum_{j:\sigma_l(j)=a}k_j,\qquad
 C_l(k)=\frac{m!}{\prod_{a\in A_l}K_{l,a}!}.
\]

For a separating set `B`, mapping a word in the `s` colors to its
layer words is injective: the tuple of symbols at each position recovers
the original color. There are `C_k` words with color counts `k` and at
most `prod_(l in B) C_l(k)` tuples of layer words with the specified
counts. Therefore

\[
                         C_k\le\prod_{l\in B} C_l(k).                  \tag{3}
\]

For every positive integer `D`, coloring `D` labelled blocks, each of
size `m`, injects into the colorings of their union. Consequently

\[
 \frac{\prod_a(DK_a)!}{(Dm)!}
       \le\left(\frac{\prod_a K_a!}{m!}\right)^D.                     \tag{4}
\]

The normalized marker contribution to a mixed moment is thus at most

\[
 \begin{aligned}
 R(k)&=\prod_l\frac{\prod_a(D_lK_{l,a})!}{(D_lm)!}\\
 &\le\prod_l C_l(k)^{-D_l}
 \le\prod_B\left(\prod_{l\in B} C_l(k)\right)^{-w_B}
 \le C_k^{-W}.                                                        \tag{5}
 \end{aligned}
\]

The middle inequality uses `C_l(k)>=1` and the capacity inequalities in
(2). This is why fractional weights are allowed. Powers of positive
rational numbers in the proof are ordinary real powers; no fractional
powers are adjoined to the coefficient field.

### 1.2. From suppression to a nonzero moment

Write `M_j=Y^(a_j)` and

\[
 H_m=\prod_i(ma_{ij})!,\qquad
 H(k)=\prod_i\left(\sum_j a_{ij}k_j\right)!.
\]

The common endpoint profile makes `H_m` independent of `j`. Integer
factorial smoothing, proved in the balanced amplification note, gives
`0<H(k)<=H_m`. Normalize the actual moments by

\[
 P_m(c)=\frac{\Lambda(F(c)^m)}{H_m\prod_l(D_lm)!}.
\]

Each pure coefficient is one, and every mixed coefficient obeys

\[
          0<[c^k]P_m=C_k\frac{H(k)}{H_m}R(k)\le C_k^{1-W}.             \tag{6}
\]

On `max_j |c_j|=1`, put `p_m(c)=sum_j c_j^m`. For `m>=2`, every
nonpure composition has `C_k>=m`, so

\[
 |P_m-p_m|\le {m+s-1\choose s-1}m^{1-W};\qquad P_1=p_1.               \tag{7}
\]

The exact comparison in Sections 1 and 3 of the balanced amplification
note applies with `W>=6s`. Its sole analytic input is the inspected
Montgomery--Vaughan formulation of Turan's second main theorem:

\[
 \max_{n\le m<n+s}|p_m(c)|
 \ge2\left(\frac{s}{8e(n+s-1)}\right)^s=:T_{s,n}.
\]

The maximum error in (7) is strictly less than
`(243/256) T_(s,n)` for every `n>=1`. Thus at least one actual moment
in the window is nonzero. The case `s=1` is immediate. Removing zero
coefficients preserves every separating set and all capacity inequalities;
`W>=6s>=6r` gives the assertion with the active count `r`. Separating
markers make the active monomials distinct even if the base monomials
coincide.

Finally, the subfield generated by finitely many coefficients over `Q`
embeds in `C`. The moments are polynomials with rational coefficients,
so a zero block over any characteristic-zero field would give the same
zero block over `C`. This proves the stated field generality.

The source and proof pins for the imported quantitative comparison are
in [its source ledger](sfc_balanced_private_amplification_source_ledger.json).
No unverified statement about the Riemann hypothesis is used.

## 2. A genuine row-merging closure

Within each layer, identify some variables by a map `pi_l:A_l -> A'_l`.
The resulting code is `sigma'_l=pi_l composed with sigma_l`. Each column
still uses exactly one variable of exponent `D_l` in that layer: variables
being merged have disjoint column supports. Suppose the final maps admit
a certificate (2).

**Corollary 2.** Every intermediate stage of these within-layer
identifications, in any order, satisfies SFC for all coefficients and
all starting orders.

Indeed, each intermediate partition refines the final partition. Every
separating set for the final code separates all intermediate codes as
well. The same weights and capacities certify every stage by Theorem 1.
The base monomials are unchanged. Different layers and base variables
are not identified with each other in this corollary.

In particular, start with a private variable for each color in every
layer. A final separating code permits the entire prescribed sequence
of identifications. This is a proved preservation statement on an
explicit class of actual integer factorial frames, rather than a
statement only about a continuous positive-measure model.

### Examples with no private variables

For four terms use the two binary layers with codewords
`00,01,10,11`. Set both exponents to `24` and give the single separating
set of two layers weight `24`. With all base monomials equal to one,
the family is

\[
 c_1X_0^{24}Y_0^{24}+c_2X_0^{24}Y_1^{24}
 +c_3X_1^{24}Y_0^{24}+c_4X_1^{24}Y_1^{24}.                            \tag{8}
\]

Each variable occurs in two terms. Every nonzero member satisfies SFC.
An arbitrary balanced monomial family can be included as the base.

Redundancy can lower the exponent per layer. Take three binary layers
with codewords `000,011,101,110`. Every pair of layers separates. Give
each of the three pairs weight `8` and every layer exponent `16`.
Each capacity is used twice, totaling `16`, while `W=24`. This gives
another four-term family with all variables shared and all-order SFC.

For arbitrary `s`, an injective binary code of length
`ceil(log_2 s)` and exponents `6s` always suffices. This construction
uses at most `2 ceil(log_2 s)` marker variables. The theorem concerns
arbitrary separating codes, not just these examples.

More generally a layer symbol may carry a monomial template with the
same positive exponent profile for every symbol in that layer. All
templates use mutually disjoint variables across layer/symbol pairs,
also disjoint from the base variables. Its capacity is the total degree
of that template: apply (4) to each coordinate and add the exponents.
Fix one common coordinate labelling of each layer's exponent profile,
and identify only equally labelled coordinates according to the prescribed
symbol coarsening; these identifications give the same corollary.
In particular, several independent unit-power coordinates can provide
the capacity instead of a single high power.

## 3. An intrinsic zero-pattern theorem

The next result is stated for the original exponent matrix itself.
It requires no construction of marker variables.

Let `A=(a_ij)` have `s>=2` columns, with identical multisets of positive
entries. Suppose its binary columns `b_ij=1_(a_ij>0)` are distinct.
Choose `ell` rows whose binary entries jointly distinguish all columns,
and put `D=min_(a_ij>0) a_ij`.

**Theorem 3.** If `D>=6s ell`, every nonzero polynomial
`sum_j c_j Y^(a_j)` satisfies SFC over every characteristic-zero field.

One can always choose `ell<=s-1`. Starting with the partition having
one block of columns, if some block contains two columns, a row
distinguishes two of them and strictly refines the partition. At most
`s-1` refinements produce singletons. Thus the simpler sufficient bound
`D>=6s(s-1)` is uniform in the number of variables and all other exponents.

### 3.1. A rowwise factorial-gap inequality

Let `g` be the convex piecewise-linear interpolation of
`g(q)=log(q!)` on the nonnegative integers. Convexity follows because
the consecutive slopes are `log(q+1)`.

Fix a composition `k` of `m` and a row `i`. Let
`K_i=sum_(a_ij>0) k_j`. The cases `K_i=0` and `K_i=m` are included in
the following inequality:

\[
 \frac1m\sum_j k_j g(ma_{ij})
       -g\left(\sum_j k_j a_{ij}\right)
 \ge\frac{K_i}{m}g(mD)-g(K_iD).                                      \tag{9}
\]

For `K_i=0`, both sides are zero. Otherwise put
`B_i=(sum_j k_j a_ij)/K_i>=D`. Conditional Jensen over the positive
entries bounds the left side below by
`h(B_i)`, where `h(x)=(K_i/m)g(mx)-g(K_i x)`.
This piecewise-linear function is continuous and, except at its finitely
many breakpoints on each compact interval, has derivative
`K_i(g'(mx)-g'(K_i x))>=0`. Hence `h(B_i)>=h(D)`, proving (9).
For `K_i=m` this last bound is zero; it does not assert equality of the
original Jensen gap when the positive entries differ.

Every column has the same number `t` of positive entries. Sum (9)
over all rows, use the common pure factorial moment `H_m`, and
exponentiate. This yields

\[
 \frac{H(k)}{H_m}
 \le \frac{\prod_i(DK_i)!}{(Dm)!^t}.                                 \tag{10}
\]

### 3.2. A distinguishing row forces suppression

Fix any row `h`. All numbers `DK_i` lie between zero and `Dm` and
their sum, excluding `h`, is `D(tm-K_h)`. Among lists of integers
bounded by `Dm` with this total, the product of factorials is at most

\[
                       ((Dm)!)^{t-1}(D(m-K_h))!.
\]

To see this, if two entries are strictly between zero and `Dm`, move
one unit from the smaller to the larger; the factorial product does
not decrease. Repeating leaves at most one interior entry and gives
the displayed expression. The endpoint cases `K_h=0,m` agree with
the same formula. Therefore (10) implies

\[
 \frac{H(k)}{H_m}
 \le {Dm\choose DK_h}^{-1}
 \le {m\choose K_h}^{-D}.                                            \tag{11}
\]

Applying the word injection to the `ell` selected binary rows gives

\[
 C_k\le\prod_h{m\choose K_h}.
\]

At least one factor is at least `C_k^(1/ell)`. Use that row in (11):

\[
                  \frac{H(k)}{H_m}\le C_k^{-D/\ell}.                 \tag{12}
\]

The normalized actual moment coefficients are consequently at most
`C_k^(1-D/ell)`. Since `D/ell>=6s`, Section 1.2 applies without any
marker factors. Active coefficient subsets retain the selected
separating rows and the same exponent bound. This proves Theorem 3.

For example, the cyclic three-term supports with columns
`(D,2D,0),(0,D,2D),(2D,0,D)` satisfy the theorem whenever `D>=36`.
All three variables are shared. The simpler triangle support
`(D,D,0),(0,D,D),(D,0,D)` is also covered.

## 4. A sharp obstruction to the near-pure estimate

The zero-pattern hypothesis reflects a limitation of the estimate,
not merely of the chosen proof. Let `a,b` be two balanced exponent
columns, and define

\[
 \beta_{a,b}=\sum_{i:a_i=0}b_i,\qquad
 C_{a,b}=\left(\prod_{i:a_i=0}b_i!\right)
                  \prod_{i:a_i>0}a_i^{b_i-a_i}>0.
\]

For the composition having `m-1` copies of `a` and one copy of `b`,
finite products of consecutive integers give

\[
 \frac{\prod_i((m-1)a_i+b_i)!}{\prod_i(ma_i)!}
       \sim C_{a,b}\,m^{-\beta_{a,b}}\quad(m\longrightarrow\infty).
                                                                         \tag{13}
\]

For `a_i>0`, the relevant factorial ratio is a finite product or
reciprocal product with leading behavior `(ma_i)^(b_i-a_i)`. For
`a_i=0` it is the constant `b_i!`. The equal profiles imply equal total
degrees, so `sum_(a_i>0)(b_i-a_i)=-beta_(a,b)`, proving (13).

Since the multinomial coefficient for this composition is `m`, a
uniform estimate `H(k)/H_m<=C_k^(-W)` with `W>0` necessarily has
`W<=beta_(a,b)` for every ordered pair. If the two zero patterns
coincide, `beta_(a,b)=0`, and no positive exponent `W` can work.
This remains true after any fixed dilation of the exponent columns.

For example, with `a=(2D,D)` and `b=(D,2D)`, the ratio in (13) is

\[
 \prod_{j=1}^{D}\frac{Dm+j}{2Dm-j+1}\longrightarrow2^{-D}>0.
                                                                         \tag{14}
\]

Every factor in (14) is greater than `1/2`, so the ratio is greater
than `2^(-D)` for every `m>=1`. Thus an arbitrarily large fixed `D`
cannot make these mixed
coefficients satisfy a power-decay estimate for all orders. This
excludes extending this particular near-power-sum proof by simply
increasing exponents. It is not a counterexample to SFC; this
example at `D=1` is already covered by the retained common-monomial
linear result.

## 5. Limits and exact verification

The binary distinction assumption is substantive. Two different positive
columns can have identical zero/nonzero patterns, and then no selected
binary rows recover their colors. Nor can arbitrary rows be treated as
partition layers: each column must select exactly one symbol in a layer.
Identifying all marker symbols destroys the separating certificate.
This may invalidate the proof while leaving SFC true for the resulting
polynomial; no SFC counterexample is asserted.

The common endpoint profile cannot be replaced merely by a common total
degree. For instance, the pure moments of `X^2` and `YZ` are `(2m)!`
and `(m!)^2`. The new results do not remove the fixed balancing factors
from the arbitrary-support enlargement in the earlier note.

[The exact checker](verify_sfc_separating_row_amplification.py) verifies
the finite certificates, literal factorial expansions, the counting and
gap inequalities on recorded finite domains, and homogeneous Macaulay
rank certificates for specified examples. It includes negative controls
for a nonseparating code, a false capacity certificate, and a failed
endpoint-profile hypothesis. Its source pins include this note and
the imported balanced amplification proof. The receipt records the
exact scope of those computations; the universal results follow from
the arguments above, not from extending the finite search range.
