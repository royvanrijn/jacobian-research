# Valuation-root models in dense residue cells

A dense additive cell has two complete-letter labels. A
bounded-degree model for powers of its exact valuation roots
in those labels forces a separation of variables. At the
generic zero of the remaining letter's sum or difference,
the two-variable divisor is reduced. This forces large
valuation orders to be fixed. The same obstruction can be
tested with a logarithmic number of grouped valuation
probes, using resultants to resolve ambiguous primes.

## RL1. Integral rational functions on dense products

Let \(\mathcal G\subseteq S_1\times\cdots\times S_r\) be a
sample of integer vectors, with \(|S_i|=H_i\), volume
\(V=\prod_i H_i\), and all coordinates of polynomial height
in \(X\). Fix the dimension and the degree bound \(D\)
before the height limit. If a reduced rational function
\(P/Q\) of degrees at most \(D\) is integral at every
point of \(\mathcal G\), then either \(Q\) is constant or
\[
              |\mathcal G|
                   \ll_{r,D,\varepsilon}
                        X^\varepsilon V/\min_i H_i .   \tag{RL1}
\]
Assume its values on the sample have a fixed polynomial
height bound, and allow its coefficients to move without
a height bound.

If \(|\mathcal G|\le O_{r,D}(V/\min_i H_i)\) there
is nothing further to prove. Otherwise evaluation is
injective on polynomials of degree at most \(2D\), by
the rectangular polynomial zero bound. That bound
follows by induction on the number of coordinates,
fixing the leading coefficient in the final variable.
Solving \(P(\mathbf t)=z(\mathbf t)Q(\mathbf t)\) on the
sample using integer minors furnishes polynomial-height
coefficients. The denominator is nonzero by injectivity
at degree \(D\). The cross product with the original
presentation vanishes identically by injectivity at
degree \(2D\), so the reduced rational function has a
polynomial-height representation. Fixed-degree integer
factor heights after reduction are polynomial by the
one-variable factor bound and a monomial-separating
substitution at fixed degrees.

If \(Q\) is nonconstant, choose a coordinate in which
its degree is positive. At all specializations of the
other coordinates off the resultant and leading
coefficient zero sets, \(P,Q\) specialize to coprime
univariate polynomials and the denominator remains
nonconstant. Their resultant in that coordinate is
a nonzero polynomial in the other coordinates since
the original polynomials are reduced. At each eligible
specialization, the nonzero value of \(Q\) at an actual
parameter divides the integer resultant by integrality
and the resultant Bézout identity. That integer has
polynomial height. Enumerating its divisors and solving
for a given nonzero value of \(Q\) costs \(X^\varepsilon\)
per specialization, with adjusted inner allowances.
On the bad specializations reinsert the whole chosen
coordinate and use the rectangular zero bound. Both
counts are at most the bound in (RL1). The argument
covers a univariate sample using a constant nonzero
resultant in place of the exceptional zero set.
[Interpolation and specialization for rational carriers (PV4)](ABC_PRIVATE_BOUNDARY_VALUATION_GRIDS.md#pv4-dense-rational-parameter-boxes-and-simultaneous-common-factors)

## RL2. Two complete labels and polynomial root powers

Put \(n=X^L\) with the height limit first at fixed \(L\).
Use genuine primitive rows \(a+b=c\), \(X/2<c\le X\),
with total original radical at most \(n^{B_0}\), and
retain any three individual caps. Use the fixed exact
valuation orders on the chart and write each letter as
a fixed full coefficient times the product of its exact
squarefree roots \(r_{T,e}^{e}\). Roots in different
slots are pairwise coprime. Fixed coefficients can
absorb any slots fixed on the chart, including a
recorded large-order tail.

Take two distinct complete roles \(U,V\). Assign their
complete values injectively to integer labels \(x,y\),
respectively. Write \(W\) for the remaining role. The
row determines an edge in \(S_x\times S_y\), and every
edge represents at most one row. Assume that, in all
remaining slots, for fixed integers \(1\le k_{T,e}\le K\),
\[
                  r_{T,e}^{k_{T,e}}=P_{T,e}(x,y),       \tag{RL2a}
\]
where the polynomials have degree at most \(D\).
Their coefficients may vary without a height bound.
Let \(S_0\) bound the sum of the degrees of the
polynomials on \(W\). A sufficient cutoff is
\[
                E_*=\max\{K,\ K\min(KD,S_0)\}.
\]
Allow any larger integer cutoff. If \(B_0 L E_*<1\)
with a fixed positive margin, the number of rows on
the chart is
\[
                      O_{L,D,K,\varepsilon}
                              ((|S_x|+|S_y|)X^\varepsilon).
                                                        \tag{RL2}
\]
All degree bounds and the slot list are fixed at the
inner height limit. The same statement holds for
reduced rational functions of numerator and denominator
degrees at most \(D\), with polynomial-height labels,
using \(S_0\) as a bound on the sum of numerator
degrees on \(W\). Apply (RL1) first to any nonconstant
denominator.

Here is the polynomial case. If the edge count is
larger than a sufficiently large degree-dependent
multiple of \(|S_x|+|S_y|\), then more than \(D\)
different \(x\)-labels have more than \(D\) neighbours,
and so do more than \(D\) different \(y\)-labels.
For those \(x\), each \(U\)-polynomial is constant as a
polynomial in \(y\), because its values are powers of a
fixed exact root of that complete \(U\)-value. Its
nonconstant-in-\(y\) coefficients vanish at more than
\(D\) values of \(x\), so all \(U\)-polynomials depend
only on \(x\). The \(V\)-polynomials depend only on \(y\).
Both lists include a nonconstant polynomial, since
the exact positive roots and fixed coefficient would
otherwise give one complete label.

Over \(\mathbb C\), adjoin variables \(z_j\) for all
slots subject to \(z_j^{k_j}=P_j(x,y)\) and impose the
original signed trinomial on those variables. Projection
to the parameter plane is finite because the power
equations are monic. Its image is closed and has
bounded degree at the fixed orders and degrees, by
Bézout and projection. If the image is proper, a
nonzero polynomial of uniformly bounded degree vanishes
on it. The two-variable grid zero bound gives
\(O_{L,D,K}(|S_x|+|S_y|)\) edges. Otherwise some component
dominates the plane and carries the signed additive
identity. Use its function field.

On that component the two complete functions \(U,V\)
belong to finite extensions of the one-variable fields
\(\mathbb C(x),\mathbb C(y)\), respectively. Each is
nonconstant: raising to a common multiple of its root
orders gives a nonconstant polynomial, as each factor
has nonnegative exponent and at least one has positive
degree. The compositum of their one-variable root fields
is the function field of the product of the smooth
normalizations of the two parameter lines in those
fields. To see this inside the component, the product is
geometrically integral, and its function field embeds
since \(x,y\) are algebraically independent. The map of
this product to the parameter plane can ramify only
above vertical or horizontal divisors. On the part
over the affine plane each complete function is regular,
because each of its roots is integral over the
corresponding affine line.

The sign choices give \(W=\pm U\pm V\) as a function.
Its zero divisor over the affine plane has no vertical
or horizontal component, by nonconstancy of the other
role. It is reduced: on any component the two projections
are dominant, so the two one-variable derivatives are
nonzero there generically. At least one nonzero partial
derivative in a local parameter prevents a repeated
component. This uses characteristic zero. In the case
\(K=1\), the same conclusion follows in the polynomial
ring because a repeated irreducible factor would divide
nonzero derivatives in the two disjoint variables.

For a nonconstant \(P_{W,e}\), take one irreducible
factor in the affine plane. Its inverse image cannot be
vertical or horizontal, because the positive exponents
on \(W\) and regularity of all slots would then force
one of those forbidden zeros. At a generic point of a
divisor over it in the product, the projection from the
product is unramified. In the larger field containing
all roots, divide the orders by the ramification index.
The resulting vanishing order of \(W\) is at least
\(e/k_{W,e}\), while reducedness bounds it by one.
Consequently
\[
                       e\le k_{W,e}
           \quad\text{for every nonconstant }P_{W,e}.  \tag{RL2b}
\]

At any divisor of the product above \(x=\infty\) with
generic \(y\), \(U\) has a pole and \(V\) is regular.
The pole of \(W\) equals that of \(U\). Comparing their
pole orders, again cancelling ramification, gives
\[
 \sum_e (e/k_{U,e})\deg P_{U,e}
          =\sum_e (e/k_{W,e})\deg_x P_{W,e}
                              \le \min(KD,S_0).
\]
There is an identical bound for \(V\) at \(y=\infty\).
The final inequality uses (RL2b) and the fact that at
most \(K\) different exact orders on \(W\) are at most
\(K\). Thus every slot of order \(e>E_*\) is a constant
polynomial on all three roles.

The positive exact roots on those slots are fixed for
the actual rows. Write \(a=A_*a_0\), \(b=B_*b_0\),
\(c=C_*c_0\) with the fixed full parts including all
orders above \(E_*\). The low parts satisfy
\[
                  a_0b_0c_0
                   \le\operatorname{rad}(abc)^{E_*}
                   \le X^{B_0 L E_*}.
\]
Bin the three sizes dyadically at height-subpower cost.
In a bin the opposite cross products of two low
pairs satisfy
\(a_0 b'_0,a'_0 b_0\le 2 a_0 b_0<C_*\) at late
heights, because \(C_* = c/c_0>X/(2c_0)\).
Primitivity makes \(A_*,B_*,a_0,b_0\) units modulo
\(C_*\). The additive relation fixes \(a_0/b_0\) there.
Rational reconstruction by the smaller cross products
then determines the exact reduced pair, hence all
three original letters. This leaves subpower many
actual rows in the identity case and proves (RL2).

In an RB4 cell of a regular balanced family, write
\(u=1/3-\alpha>0\). The two complete dictionaries have
size \(n^{u+o(1)}\) and the cell has \(n^{2u-o(1)}\)
rows. On that cell, a covering by at most \(n^{b}\)
two-label models as above, including the price of fixing
their tails, is impossible whenever \(b<u\) with fixed
slack. Sum (RL2) over the charts using the original
two dictionary sizes of the cell. No sum over the
other residue cells is used. This includes all fixed
positive \(u\) in the balanced dictionary range.
At \(u=0\) the two dictionary sizes give no power gain.
[Balanced residue cells (RB4)](ABC_BALANCED_DICTIONARY_RESIDUE_CELLS.md#the-actual-residue-cell-sizes)

## RL3. Separate modular rational models

The finite-minor lift can produce (RL2a), or its
rational analogue, from modular descriptions of
powered roots or other integer labels. Let \(\mathbf t\)
be \(r\) common short integer parameters of magnitude
at most \(n^v\) after translation. Let \(z\le n^C\)
be an integer label. On the chart assume
\[
                  P(\mathbf t)\equiv z Q(\mathbf t)
                                \pmod M,\qquad M\ge X^\nu ,
\]
with \(Q(\mathbf t)\) a unit on every actual row and
degrees at most \(D\) over \(\mathbb Z/M\mathbb Z\).
Let \(s\) bound the number of monomials of degree at
most \(D\). If
\[
                          Ls(C+3Dv)<\nu
                                                        \tag{RL3}
\]
with a fixed positive margin, a bounded refinement at
the fixed \(L,r,D\) makes \(z\) an exact rational
function of the common parameters of degrees at most
\(D\). Applying the argument to a finite list of labels
permits different moduli, with no extra normalized
price at fixed \(L\). The labels can be \(z=r_{T,e}^k\),
using their powered numerical height in \(C\).

Choose a rational basis \(B\) of all evaluated monomial
columns and a nonzero square minor \(\Delta_0\) of this
basis. Write \(s_0\le s\) for its size. Multiplying each
monomial column by \(\Delta_0\) places it in the integral
span of \(B\). Consequently the modular relation gives
an integral relation on \([B,zB]\) modulo \(M\), whose
denominator part evaluates to \(\Delta_0 Q(\mathbf t)\).
The greatest common divisor of all its coefficients
with \(M\) divides \(\Delta_0\), since \(Q\) is a unit
on the sample. Multiply the relation by the adjugate
of any full minor \(\Delta\) of \([B,zB]\). This yields
\(M\mid \Delta\Delta_0\). The product has absolute
value at most \((2s_0)!s_0! n^{s_0(C+3Dv)}\), so every
full minor vanishes.

An exact rational relation from the evaluation kernel
now has nonzero denominator as an evaluated function,
since the \(B\) columns are independent. Coefficients
can be taken at polynomial height using minors.
On the rows where it vanishes, the monomial evaluation
rank drops. Reapply the argument there with a smaller
basis. After at most \(s\) rank drops the refined rows
all have exact rational models on the part where their
new denominator is nonzero. All new equations use the
same original common parameters. Moduli and auxiliary
coefficient lifts need no common prime support.
[Polynomial modular lift (PV5)](ABC_PRIVATE_BOUNDARY_VALUATION_GRIDS.md#pv5-exact-rational-carriers-from-separate-deep-moduli)

## RL4. A logarithmic family of grouped probes

Keep all the exact valuation roots not fixed on the
chart, indexing them by \(j=1,\ldots,m\). For a fixed
matrix \(A=(A_{\ell j})\) of nonnegative integers
use probe labels
\[
                         M_\ell=\prod_j r_j^{A_{\ell j}} .
\]
Assume that the columns are nonzero and lie on distinct
rational rays. Write \(A_{\bullet j}=g_j u_j\) with
\(u_j\) primitive and put \(K=\max_j g_j\).
These probe labels determine all the exact roots on
a row, by the valuation vector of each prime. Suppose
the probes are given by polynomial or reduced rational
functions of common integer parameters in a product
box of volume at most \(n^v\) and shortest side at
least \(n^\kappa\). Degrees and the number of
polynomials are fixed at each \(L\). Let \(S_*\) bound
the sum of their numerator degrees. One may use
polynomial-height coefficients, or obtain them by
the interpolation argument of (RL1) off a chart
with \(O(n^{v-\kappa})\) points. This requires
polynomial height for the probes, as holds for every
matrix fixed at \(L\) and original height at most \(X\).

For any fixed \(\chi>0\), the chart has a part of size
at most \(n^{v-\min(\chi,\kappa)}X^\varepsilon\), and
the remaining rows are covered by \(n^\chi X^\varepsilon\)
refinements with polynomial power models
\[
                              r_j^{g_j}=c_j G_j(\mathbf t)
                                                        \tag{RL4}
\]
for fixed nonzero rational \(c_j\), where the sum of
the degrees of the integer polynomials \(G_j\) is
at most \(S_*\). Different refinements change the
constants. The proof uses the actual squarefree
and pairwise coprime exact valuation roots, including
coprimality within an original letter.

Use the primitive irreducible factors across all the
probe numerators and denominators. On an actual sample
apply the common shear, resultant and CRT count of
PV4 to their simultaneous collision primes. Exclude
contents and the leading coefficients in the shear
from that product, treating them as a fixed prime pool.
Parameters for which the collision radical exceeds
\(n^\chi\), or the specialized resultant vanishes, give
the first asserted count. On the rest, fix the
collision support at \(n^\chi\) cost and add it to the
pool. At each prime outside the pool, at most one
of the primitive factors \(F\) can vanish at the
parameters. A vanishing factor cannot occur in
any denominator of the reduced presentations, by
integrality of the integer probes. If \(w_F\) is
its vector of multiplicities in all numerators,
then the valuation vector of the probe values at
that prime is exactly \(v_p(F(\mathbf t))w_F\).
If this vector is nonzero it equals \(A_{\bullet j}\)
for a unique \(j\), since \(r_j\) is squarefree and
the slots are pairwise coprime.

Assign \(F\) to slot \(j\) if \(w_F=h_F u_j\) for a
positive integer \(h_F\) and \(F\) lies on no
denominator. Let \(G_j\) be the product of the assigned
\(F^{h_F}\). Their values off the pool have exactly
the prime powers of \(r_j^{g_j}\) off the pool.
Factors that have no assignment have all their
actual prime support on the pool. A given factor
is assigned to at most one \(j\), and the sum of
the \(h_F\deg F\) is bounded by the sum of all
numerator degrees. On each prescribed pool fix
the supported parts of both sides at height-subpower
cost, as they have polynomial height at fixed \(L\).
Fix their signs as well. This gives (RL4).
[Simultaneous collision-prime count (PV4)](ABC_PRIVATE_BOUNDARY_VALUATION_GRIDS.md#pv4-dense-rational-parameter-boxes-and-simultaneous-common-factors)

One may choose \(K=1\) and only
\(\lceil\log_2(m+1)\rceil\) probes, using distinct
nonzero binary signatures for the slots. Separating
the roles first gives the explicit probes on each
letter
\[
                M_{T,i}=\prod_j
                        r_{T,j}^{\lfloor j/2^i\rfloor
                                       -2\lfloor j/2^{i+1}\rfloor} .
\]
With cutoff \(B/L\) these use \(O(\log(1/L))\) probes
per role. They are the ratios of consecutive
power-divisor roots \(F_{T,i}/F_{T,i+1}^2\), where
\(F_{T,i}=\prod_p p^{\lfloor v_p(T)/2^i\rfloor}\),
after excluding the fixed higher-order tail.
These operations require no individual layer
interpolation. Distinct rays matter for the
factor argument: with two parallel columns \(1,2\),
the mask \(t\) can receive a prime from the first
slot when \(t=p\) and from the second when \(t=p^2\).

For two label parameters as in RL2, and
\(B_0 LK\max(1,S_*)<1\) with slack, sum (RL2)
over the retained probe refinements. The per-chart
count is at most
\[
       X^\varepsilon\bigl(n^{v-\min(\chi,\kappa)}
                               +n^\chi(|S_x|+|S_y|)\bigr).
                                                        \tag{RL4b}
\]
In particular in a balanced cell of positive exponent
\(u\), a box of side at most \(n^{u+o(1)}\) in both
labels gives a strict saving for \(0<\chi<u\)
by fixed margins, including enough room for tail
and model catalogue prices. Probe degrees may
be bounded while their number grows logarithmically,
since then \(LS_*\to0\) for binary probes.
The modular criterion RL3 applies to each probe at
its own numerical modulus.

## RL5. Effective Hilbert irreducibility in a two-label cell

It suffices for a weaker but strict saving to have one
nonzero bounded-degree polynomial relation for each exact
root, without a specified rational function or a chosen
branch. Use the same two label parameters, now of size at
most \(H\) after translation, and assume
\[
                         F_j(x,y,r_j)=0
\]
in each slot for nonzero \(F_j\in\mathbb Q[x,y,Z]\)
of bounded total degrees and parameter degrees at most
\(D\), with arbitrary coefficients at the fixed degrees.
Here parameter degree is the total degree in \(x,y\)
over the coefficient ring \(\mathbb Q[Z]\). If
\(B_0 LD<1\) by a
fixed positive margin, the count in this chart is at
most \(H^{3/2} X^\varepsilon\), up to constants at
the fixed \(L\), degree bounds and slot list. Here \(H\ge2\)
is polynomial in \(X\), and the roots are the exact
integral valuation roots as before.

First choose one relation for each slot from the actual integer monomial evaluation kernel
using the original degree bounds in the root and the
parameters. The
parameters and roots have polynomial height, so these
relations have polynomial-height coefficients by integer
minors. Factor them over \(\mathbb Q\). Factors with zero
root degree, and specializations on which the leading
root coefficient vanishes, give \(O_L(H)\) parameter pairs
per factor by the grid bound. An irreducible factor of root
degree at least two that specializes to have an integer root
at any other parameter pair specializes reducibly.
Theorem 1.10 of Cluckers, Dèbes, Hendel, Nguyen and
Vermeulen, with two parameters and one root variable,
bounds these specializations by \(O_L(H^{3/2}X^\varepsilon)\),
with dependence on the fixed total degree bounds.
Its polynomial dependence on logarithmic coefficient
height is absorbed by \(X^\varepsilon\), since
Kronecker substitution also bounds the factor heights
at fixed degrees. Sum over the finitely many factors
and slots at the fixed budget.

On the remaining pairs select one root-linear factor
in each slot, from a finite list. The parameter degree
cannot increase upon taking a factor, by additivity of
total degree in \(x,y\) over \(\mathbb Q[Z]\). Thus
these give rational root models of degree at most \(D\), and the leading coefficients
are nonzero on the chosen parameters. Apply RL1 and
RL2 at \(K=1\) to count at \(O_D(HX^\varepsilon)\).
At an RB4 cell with \(H\le n^{u+o(1)}\) and \(u>0\),
the joint model and tail prices can therefore be any
fixed amount below \(u/2\) with a margin. A separate
deep congruence in each root and the common label
parameters supplies the exact relations if each
coefficient vector is primitive modulo its modulus
and the full augmented monomial determinants are
smaller than that modulus. Indeed multiplying any
full monomial square matrix by the modular coefficient
vector and applying the adjugate makes its determinant
zero modulo the modulus, and then zero over the
integers. The resulting evaluation kernel supplies one
exact relation per slot. The moduli can vary by slot.
[Cluckers–Dèbes–Hendel–Nguyen–Vermeulen, Theorem 1.10](https://arxiv.org/pdf/2311.16871)

A useful version lifts a numerical powered factor
\(z=r_j^{k_j}\le X^{h_j}\) from a polynomial
congruence of degree \(D\) on short parameters.
With a basis of \(s_0\) degree-\(D\) monomial columns
and parameter magnitudes at most \(n^v\), use just
one appended \(z\) column. Its augmented minors
vanish if
\(h_j+Ls_0Dv<\log_X M_j\) by a fixed margin.
The lifted polynomial has parameter degree at
most \(D\), so \(r_j^{k_j}-P_j(\mathbf t)=0\) can
enter RL5 even at large \(k_j\) fixed at the inner
height limit. For a modular rational function with
unit denominator, use the larger determinant height
in RL3. A missing slot may also be supplied by a
fixed exact value, charging the joint catalogue of
all such records once.

The same Hilbert step can be applied to relations for
the grouped probes of RL4. On the root-linear
specializations use the rational models of the probes
and the collision refinement there. In a balanced
cell of side at most \(n^{u+o(1)}\), take
\(\chi=u/2\). The Hilbert-exception branch, the
large-collision branch, and the polynomial-power
branch all cost at most \(n^{3u/2+o(1)}X^\varepsilon\),
provided \(B_0 LK\max(1,S_*)<1\) by a fixed margin.
For binary probes the number of separate exact
algebraic relations needed is logarithmic in the
number of active slots. Both criteria keep the
third original radical cap and allow arbitrary
original least-repair eligibility.

## The missing information in complete registers

The preceding lifts require a relation for a root or
a grouped probe. Complete register equations by
themselves have a free torus of hidden directions:
for two slots of orders \(e,f\) on one role, the
transformation
\[
 (r_e,r_f)\longmapsto
        (r_e t^{f/g},r_f t^{-e/g}),\qquad g=\gcd(e,f),
\]
preserves that entire complete value, hence all native
complete registers. Over an algebraic closure the
monomial map from \(m\) nonunit slots to the three
complete role values has geometric generic fibre
dimension \(m-b\), where \(b\le3\) is the number of
roles with moving slots. The transformation also acts
on units modulo any common numerical modulus. Thus
elimination of the complete register equations
with no extra arithmetic input leaves those directions
free in the algebraic model. The transformations do
not preserve positive squarefree pairwise-coprime
integer roots at the given heights. Extracting the
separate probe relations on the genuine rows still
requires that arithmetic information. At occurrence
\(u=0\), a chart spanning multiple residue cells
or a different mechanism is also needed.

Walsh's inverse sieve gives a bounded-degree relation
from sufficiently small residue occupancy of an integer
set across many primes. A fixed complete-register cell
provides congruences on the three complete values, so
that theorem's occupancy hypothesis would have to be
verified for the joint label and probe coordinates
themselves before using it for RL5.
[Walsh, Theorem 1.1 and Theorem 2.4](https://arxiv.org/pdf/1105.1551)

## RL6. A squarefree anchor and a complete-letter model

One extra probe removes the collision price entirely.
Include the total retained radical \(R=\prod_j r_j\)
as a modeled integer label. Suppose that the other
probe columns of \(A\) are distinct, including a
possible zero column, and fix \(A_\ell\ge\max_j A_{\ell j}\).
There is no distinct-ray assumption in this version.
On a two-label chart with \(H_1,H_2\) actual labels,
either the row count is \(O_L((H_1+H_2)X^\varepsilon)\)
or every probe model can be taken polynomial by RL1,
also when originally given by rational functions.
Denote those polynomials by \(P_R,P_\ell\). Since
\(R^{A_\ell}/M_\ell\) is integral on the sample,
apply RL1 again to \(P_R^{A_\ell}/P_\ell\).
Its values are of polynomial height at each fixed
\(L\). Unless this gives the same small count,
\(P_\ell\) divides \(P_R^{A_\ell}\) over \(\mathbb Q[x,y]\).

Enlarge a fixed pool by the contents and denominators
in the primitive factorizations. For a prime outside
it, when \(p\mid R\) the squarefree value of \(P_R\)
means that exactly one primitive irreducible factor
\(F\) can contribute, necessarily with polynomial
multiplicity one and \(v_p(F(x,y))=1\). Its
multiplicity vector across the \(P_\ell\) specifies
the unique slot \(j\), by distinctness of the probe
columns. All factors with an invalid vector or with
multiplicity greater than one in \(P_R\) have all
their actual prime support on the fixed pool.
The others give polynomial models
\(r_j=c_jG_j(x,y)\) after fixing the supported parts
and signs at height-subpower cost. Their total
degree is at most \(\deg P_R\). The primitive
root powers need no collision catalogue, since
the squarefree total already forbids a collision
outside the pool. Apply RL2 with \(K=1\) if
\(B_0 L\max(1,\deg P_R)<1\) with a fixed margin.
This proves the small count on the two-label chart.

For \(m\) slots use \(R\) and two digit probes by
encoding distinct slot numbers with two digits in
base \(B=\lfloor\sqrt m\rfloor+1\), increasing \(B\)
by one if necessary. Each digit probe has numerical
height at most \(n^{(B-1)B_0}\). With truncation
\(m=O(1/L)\) it is \(X^{O(\sqrt L)}\). Thus three
integer labels suffice as input to RL3 on short
parameters with fixed positive modulus depths.
Alternatively binary digits give \(O(\log(1/L))\)
labels of height \(X^{O(L)}\). The already known
least multiplier on a single retained role whose
orders are at most \(M\) is another valid tag:
\(g_M(T)=\prod_e r_{T,e}^{M-e}\). Together with that
role's radical it gives all its exact orders by the
same argument, with the multiplier's actual height
used in a modular lift.

There is a two-model version. Suppose the total
retained radical and the remaining complete letter
\(W\), or that complete letter divided by its known
full coefficient, are given by rational functions
of the two label parameters. Assume their degrees
are bounded at the fixed \(L\), the labels have
polynomial height, and \(B_0 L\max(1,\deg P_R)<1\)
with slack for the reduced numerator of the radical
model. Then the two-label chart has at most
\[
                   O_L\bigl(((H_1H_2)^{2/3}+H_1+H_2)
                                      X^\varepsilon\bigr).       \tag{RL6}
\]
To prove this, use RL1 to obtain a polynomial \(P_W\)
for \(W\), or for the known normalization and then
multiply by its fixed coefficient. If \(P_W\) is not a sum of an \(x\)-polynomial
and a \(y\)-polynomial, apply the point–curve bound
with two degrees of freedom. Absorb the original
signs into \(U,V\). Use the graphs
\(w=P_W(x,y)-U(x)\) indexed by the \(x\)-labels,
and the points \((y,V(y))\). Some coefficient of
a positive power of \(y\) is nonconstant in \(x\),
so each graph repeats at most \(O_L(1)\) times.
Distinct graphs meet at most \(\deg P_W\) times.
Some coefficient of a positive power of \(x\) is
also nonconstant in \(y\). Partition the \(y\)-labels
into at most \(\deg P_W\) lists on each of which
that coefficient takes distinct values. For two
points in such a list, \(P_W(X,y_1)-P_W(X,y_2)\)
is a nonconstant polynomial in \(X\). Thus the two
points belong to at most \(\deg P_W\) graphs by
subtracting their incidences. Pach–Sharir applies
to each list and gives (RL6).
[Pach–Sharir incidence bound, stated with its two degrees of freedom](https://math.nyu.edu/~pach/publications/DistancesOnCurves101113.pdf)

If \(P_W\) is separated, the additive relation makes
the two complete labels equal to their assigned
one-variable polynomials plus opposite fixed offsets
on each connected component of the edge graph,
with the corresponding fixed signs. The three
normalized complete values \(T/T_*\) are therefore
polynomial probes on that component. Their tag
columns identify the original role and exact order.
They have bounded degree at each \(L\), and their
quotients from a sufficient power of \(R\) are
integers on the actual rows. The squarefree anchor
argument above gives the linear-size count in the
two label counts of that component. Components have
disjoint label sets on both sides, so there is no
additional catalogue cost when summing them.

In particular, use the original complete \(U,V\)
as the parameters themselves. Then \(W\) is
already a separated linear function and the only
hypothesis is a rational model for \(R\) in these
two complete values. In this case the chart has
\(O_L((H_1+H_2)X^\varepsilon)\) rows. This version
does not need short numerical parameters, since the
divisibility and separation arguments depend only
on the label counts and polynomial heights.

## RL7. A shared index across the occurrence range

The following is a sufficient extraction that also
covers \(u=0\). On balanced charts with the RB4
exponents, write \(u=1/3-\alpha\) and \(s=1/3-2u\).
Suppose one can index the retained cells by
\(z\in[1,n^{s+o(1)}]\) and inside each cell index the
two free labels by \(x,y\in[1,n^{u+o(1)}]\). The label
indexing may vary with the cell. Suppose the complete
exact-root tuple is given by bounded-degree rational
functions in the common three parameters on each
model, after fixing an exact high tail. Let \(b\)
bound the combined price of the model records and
tails, excluding the three index coordinates.
The bounds on degrees are fixed before the budget
decreases. Choose a compact positive-order diffuse
band of mass at least \(\rho\) with maximum slot
price tending to zero, as in PV4.

If \(u>b\) with sufficient slack, fix the cell
index and use RL2 with \(K=1\). This costs at most
\(n^{s+u+b+o(1)}\) rows in total, with arbitrarily
small normalized allowances. Alternatively fix
the two free-label indices across cells, at price
\(2u+o(1)\). The surviving model has just one
short integer parameter. A nonconstant denominator
gives height-subpower count by RL1. For polynomial
functions apply PV6 to the original roots and their
three original caps. This gives at most
\(n^{(1-\rho)/3+b+2u+o(1)}\) rows in the balanced
limit. Exact coefficient and eligibility records
may be kept in both arguments. Thus \(b<\rho/9\)
by a fixed margin makes at least one route strictly
smaller than \(n^{1/3}\), uniformly in
\(u\in[0,1/6]\). Use \(u>b+\eta\) or
\(u\le b+\eta\) with one common \(\eta\) small enough
for both margins. GE1 allows \(\rho\) arbitrarily
close to \(1/2\), giving the catalogue threshold
\(b<1/18\) with fixed slack. This argument requires
one genuine common index model over the cells.
[Polynomial short-parameter exclusion (PV6)](ABC_PRIVATE_BOUNDARY_VALUATION_GRIDS.md#pv6-short-polynomial-parameters-of-arbitrary-box-volume)
