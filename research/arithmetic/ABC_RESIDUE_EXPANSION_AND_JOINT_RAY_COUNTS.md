# Multiplicative expansion in absolute residue cells and joint ray counts

A deep absolute modulus can charge a list of residues by its
relative multiplicative expansion. This strengthens the
square-root cardinality price used for balanced abc cells.
In particular, every sufficiently large part of the occupied
residue list at a positive dictionary surplus has nearly
maximal growth under every bounded number of modular products.
The same Fourier count applies to a joint ray group for several
actual primitive rows. A common numerical target divisor can
then charge shared target primes at half their merged radical
price.

## The relative expansion price for a deep absolute modulus

**RX1.** Put \(n=X^L\). Fix \(B,w>0\) and \(r,b\ge0\).
Suppose \(X^w\le D\le X^B\), and write \(G=(\mathbb Z/D\mathbb
Z)^\times\). Let \(S,\Theta\) be nonempty subsets of \(G\),
with \(|S\Theta|\le n^b\). For every \(\eta>0\) and sufficiently
small fixed \(L>0\), the number of positive integers \(V\le X^B\) with
\(\operatorname{rad}(V)\le n^r\) and \(V\bmod D\in S\) is
\[
       O_{B,w,r,b,\eta,L,\varepsilon}\left(
          n^{r/2+\eta}\left(\frac{|S\Theta|}{|\Theta|}\right)^{1/2}
          X^\varepsilon\right)
                     \qquad(\varepsilon>0).
\]
A prescribed polynomial-height prime support is allowed at
zero price in \(r\), with the constant also depending on a
fixed height exponent of the support. The auxiliary set
\(\Theta\) can depend on the chart and on \(S\). Taking a
singleton gives the square-root cardinality estimate.

Here is the Fourier step independently of the modulus.
On any finite abelian group with counting convolution, use
the nonnegative majorant of \(\mathbf1_S\)
\[
                F=\frac{\mathbf1_{S\Theta}*
                          \mathbf1_{\Theta^{-1}}}{|\Theta|}.
\]
For each \(s\in S\), every \(\theta\in\Theta\) contributes the
term \(s=(s\theta)\theta^{-1}\). Take the Fourier transform as
the unnormalized sum against a character, and use probability
counting measure on the dual. Parseval and Cauchy–Schwarz give
\[
    \|\widehat F\|_1\le
         \left(\frac{|S\Theta|}{|\Theta|}\right)^{1/2},
       \qquad  \|\widehat F\|_\infty\le |S\Theta|.
\]
For \(0<\xi\le1/2\), interpolate at \(p=(1-\xi)^{-1}\), with
conjugate exponent \(1/\xi\). The relative expansion is
at least one, so
\[
  \|\widehat F\|_p\le
           \left(\frac{|S\Theta|}{|\Theta|}\right)^{1/2}
            n^{b\xi}.
\]
For character evaluations after a homomorphism with
kernel of size \(K\), the \(L^\infty\) bound is unchanged
and the \(L^1\) bound acquires a factor at most \(K\).
Indeed the squared \(L^2\) norm of the pushforward
of any \(\mathbf1_A\), counting multiplicities, is at
most \(K|A|\). This covers bounded and height-subpower
kernel sizes.

To complete the integer count, first fix the powers on the
prescribed support at height-subpower cost. Fix the remaining
full part on valuations above \(J\) at cost
\(X^{B/(J+1)+o_X(1)}\), choosing \(J\) at the fixed budget so
that \(B/(J+1)<\eta L/8\). On the rest write \(V=d\prod_{e\le
J}v_e^e\), with squarefree \(v_e\) and with \(d\) fixed and a
unit. Bin the sizes \(v_e\le n^{\theta_e}\) in a fixed finite
mesh, with \(H=\sum_e\theta_e\le r+\eta/4\), and omit zero
slots. The target list is \(d^{-1}S\), a translation that
preserves its expansion. For any positive slot use
\[
       m_e=\left\lceil\frac{H}{2\xi\theta_e}\right\rceil,\qquad
       F_e(\chi)=\sum_{v_e}\chi(v_e^e).
\]
Then \(\sum_e(2m_e)^{-1}\le\xi\) and
\(m_e\theta_e\le(1+1/(2\xi))H\).
In the \(2m_e\)-th character moment the two products of
\(m_e\) roots have equal \(e\)-th powers modulo \(D\). Fix the
first product. The second lies in at most \(X^{o_X(1)}\)
unit residue classes, by the \(e\)-torsion bound in
\((\mathbb Z/D\mathbb Z)^\times\). For fixed \(e\), a prime power
contributes at most \(2e\) roots even when it divides \(e\).
Since \(L(1+1/(2\xi))H<w\) at sufficiently small fixed \(L\),
each second product is below \(D\), so each such class
determines one exact product. The ordinary divisor estimate
bounds its ordered factorizations by \(X^{o_X(1)}\) at fixed
\(m_e\). It follows that
\(\|F_e\|_{2m_e}\le n^{\theta_e/2}X^{o_X(1)}\).
Fourier inversion, the majorant and Hölder give the product
of these bounds times \(\|\widehat F\|_p\). Choose \(\xi>0\)
first, in terms of \(b,\eta\), to make \(b\xi<\eta/4\).
The tail, mesh and height-subpower counts fit the remaining
allowances. The case \(H=0\) needs only the single fixed
integer after the initial enumeration. No threshold on \(L\)
depends on the eventual \(J\), because the product height
uses \(H,\xi\) and the root-of-unity cost is subpower at each
fixed \(J\).
[The cardinality-price estimate (RB2)](ABC_BALANCED_DICTIONARY_RESIDUE_CELLS.md#regular-balanced-cohorts-and-absolute-residue-pricing)
· [Fixed-support counts (EB3)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

## Every large balanced residue list has maximal finite-order growth

**RX2.** Use a regular balanced cohort of exponent \(1/3-o(1)\)
and dictionary exponent \(\alpha\in[1/6,1/3]\), in the
height-first convention of RB1–RB4. Write
\(u=1/3-\alpha\) and \(s=2\alpha-1/3\). Fix a cheap full
divisor \(D\) on an eligible original letter \(W\) and make the
RB4 deletions on the resulting chart. Given any positive
tolerance, the price of the chart and the RB4 tolerances may
be taken sufficiently small before going to sufficiently
late fixed budgets. Write \(S\) for the set of retained
absolute residues of the free letter \(U\). For every fixed
bound \(b\), any \(T\subset S\) and any nonempty auxiliary
\(\Theta\) with \(|T\Theta|\le n^b\) satisfy
\[
    \frac{|T\Theta|}{|\Theta|}
         \ \ge\ n^{-s-o(1)} |T|^2 .                         \tag{RX2}
\]
The \(o(1)\) is in the stated outer sense and uniform in
the two sets once their height and size bounds are fixed.
The same statement holds for \(S^{-1}\). More generally,
any subcollection of residues retaining chart rows of
exponent \(1/3-o(1)\) has expansion at least \(n^{s-o(1)}\)
against every admissible \(\Theta\), even without the RB4
lower bound on each individual cell.

Indeed each retained cell witnesses \(n^{u-o(1)}\)
distinct original \(U\)-values. They are disjoint over
the residues \(T\), hence RX1 and the original
\(n^{1/3+o(1)}\) cap give
\[
         |T| n^{u-o(1)}
              \ \le\ n^{1/6+o(1)}
                          (|T\Theta|/|\Theta|)^{1/2}.
\]
Now \(2u-1/3=-s\). Inversion replaces an auxiliary set
by its inverse. For the variant without per-cell lower
bounds, divide the number of retained rows by the
original maximum degree \(n^{u+o(1)}\) on \(U\) and apply
RX1 to its at least \(n^{\alpha-o(1)}\) original values.
This proof uses the original three-cap cohort through
RB1–RB4. The price calculation itself uses the complete
letter cap and the full divisor, with no valuation
profile assumed for the free letter.

In particular if \(|T|=n^{s-o(1)}\), then all fixed
mixed product sets \(T^a T^{-b'}\), for nonnegative
integers \(a,b'\) with \(a+b'\ge1\), have size
\(n^{(a+b')s+o(1)}\). Start with a singleton and
apply the expansion inequality to each multiplication,
using \(T^{-1}\) when required. The successive sets
all have polynomial size in \(n\) since the order
is fixed before decreasing the budget. The opposite
upper bound follows from \(|S|\le n^{s+o(1)}\).
It also applies to different subsets of
size \(n^{s-o(1)}\) at the successive factors, and
to replacing any factor \(T\) by its pointwise
power set \(T^{[e]}=\{t^e:t\in T\}\), for fixed
nonzero \(e\). For positive \(e\), partition an
auxiliary \(\Theta\) into cosets of \(G^e\). In
each coset write the elements as \(z\widetilde
\theta^{\,e}\) for a fixed coset representative
\(z\) and one selected \(e\)-th root per element.
On this set the pointwise power map can shrink the
product \(T\widetilde\Theta\) by at most
\(|G[e]|=X^{o_X(1)}\). This is the same uniform
unit torsion estimate as in RX1. If the image
product is within the prescribed bound \(n^b\),
the lifted product is within \(n^{b+1}\) at all
sufficiently large heights at fixed \(L\).
Apply RX2 in each lifted coset and sum over the
disjoint image cosets. Negative \(e\) also use
inversion. Thus the growth exponent \(ms\)
holds for every product of \(m\) full-size
subsets with arbitrary prescribed nonzero
integer powers at its slots.

This gives a saturation obstruction:
a power-large loss in growth at any fixed order
in such a chart rules out that regular balanced
cohort, provided the loss exceeds the actual chart
and pruning tolerances. For the minimum dictionary
\(s=0\), the expansion bound has no positive
normalized exponent. For any fixed
\(\alpha>1/6\), it excludes low-doubling residue
structures in full-size lists, including occurrence
\(\alpha=1/3\).
[Actual residue-cell sizes (RB4)](ABC_BALANCED_DICTIONARY_RESIDUE_CELLS.md#the-actual-residue-cell-sizes)

## Signed residue words and diagonal energy on the cells

**RX7.** Fix \(m\ge2\), nonzero integer weights
\(e_1,\ldots,e_m\), heights and a common deep modulus
as in RX1. Count tuples of positive integers \(V_i\)
with individual radical caps \(n^{r_i}\), each a unit
modulo \(D\), for which \(\prod_i V_i^{e_i}\bmod D\)
lies in a list \(S_0\) of size \(R\).
For every fixed positive allowance the count is
\[
          O\bigl(n^{(\sum_i r_i)/2+\eta}
                         R^{1/2} X^\varepsilon\bigr).
\]
If \(|S_0\Theta|\) is bounded by a fixed power of \(n\),
the factor \(R^{1/2}\) can instead be
\((|S_0\Theta|/|\Theta|)^{1/2}\).
All original integers are counted, including
repetitions. A prescribed polynomial-height prime
pool on the tuples can be removed from the radical
caps.

For the proof, fix the powers on the pool and
the entire part supported on primes of any
valuation above an inner cutoff \(J\), at
subpower and small tail costs as in RX1.
The height bound for the product of all columns
is \(X^{mB}\), so the union radical of that
tail is at most \(X^{mB/(J+1)}\).
Group remaining primes by joint exact
valuation vectors \(v=(v_1,\ldots,v_m)\).
A nonzero vector has a squarefree root record
and word exponent \(\lambda(v)=\sum_i e_i v_i\).
Let \(b_v\) be the normalized radical mass in
that vector. If \(\lambda(v)=0\), at least two
coordinates are positive. Fix the complete
zero-exponent records at full price
\(\sum_{\lambda(v)=0}b_v\). All other slots
have nonzero, possibly negative, word exponents.
The amplified finite-group moments from RX1
apply to them. For a negative exponent,
equality of powers of unit products is
equivalent to equality after inversion.
Charge half the remaining sum
\(\sum_{\lambda(v)\ne0}b_v\), then the chosen
Fourier price for the list. Use a fixed finite
mesh at each budget. The exact slot records
determine the tuple. A prime in a zero slot
has been paid at no more than half its
original incidence, and this is also true at
every other prime. More precisely, this saves
half the mass of repeated incidences beyond
one in the nonzero slots and beyond two in the
zero slots, when these masses are stabilized
before taking the height limit.

In the setting of RX2, each retained residue
has at least \(n^{u-o(1)}\) original complete
\(U\)-values. Lift each of the \(m\) residue
slots independently to these values, preserving
its residue and the original cap \(1/3+o(1)\).
Thus the number of residue words with output in
any specified list of size \(R\) is at most
\[
                           n^{ms/2+o(1)}R^{1/2},
                  \qquad s=2\alpha-\tfrac13 .
\]
For example, the positive product moment counting
ordered pairs of \(j\)-tuples with the same product
is \(n^{js+o(1)}\) on a full-size retained list.
The upper bound uses the single output \(1\) for
the word with \(j\) plus and \(j\) minus signs.
The lower bound uses identical tuples. This is
diagonal at the exponent scale at every fixed
order. More generally, any selected collection of
\(m\)-tuples of residues of size
\(n^{ms-\rho-o(1)}\) has signed image of
size at least \(n^{ms-2\rho-o(1)}\).
The moment order remains fixed before the
small-budget threshold.
[Complete-letter witnesses in residue cells (RB4)](ABC_BALANCED_DICTIONARY_RESIDUE_CELLS.md#the-actual-residue-cell-sizes)

## Joint target ray groups and shared radical prices

**RX3.** Fix an integer \(m\ge2\), \(w,B>0\), and \(g,\ell,b\ge0\).
Use \(m\) genuine positive primitive rows \(a_i+b_i=c_i\le X\)
with all three original radical caps imposed. Prescribe
\(K_i=\mathbb Q(\sqrt{-D_i})\), where
\(D_i=\operatorname{sf}(a_i b_i)\) is the squarefree kernel,
and suppose \(D_i\le n^g\). At least two fields are distinct.
Prescribe an odd \(M\ge X^w\) dividing every \(c_i\), and a
positive integer pool \(P\le X^B\).
Put \(w_i=\prod_{p\nmid M}p^{v_p(c_i)}\) and prescribe a
list \(S\) of the vectors of normalized residues
\[
             (a_i/w_i\bmod M)_{i=1}^m
          \quad\text{in }\
             U=((\mathbb Z/M\mathbb Z)^\times)^m/
                      \operatorname{diag}((\mathbb Z/M\mathbb Z)^\times).
\]
For \(p\nmid PM\), let \(s(p)\) be the number of distinct
positive values in \((v_p(c_i))_i\). Assume
\[
                        \prod_{p\nmid PM}p^{s(p)}\le n^\ell .
\]
Then for every \(\eta>0\) and sufficiently small fixed \(L\),
the ordered tuple count is at most
\[
       O\left(n^{\ell/2+\eta}
               \min\left\{|S|,\frac{|S\Theta|}{|\Theta|}\right\}^{1/2}
                  X^\varepsilon\right)
\]
for any nonempty \(\Theta\subset U\) with \(|S\Theta|\le n^b\).
The cardinality arm holds without a size bound on \(S\).
All constants are uniform in the fields and prescribed data
at the given caps. No radical bound on a new shift is needed.
For two rows the charged quantity is the product of the
outside radicals of \(\gcd(c_1,c_2)\),
\(c_1/\gcd(c_1,c_2)\), and \(c_2/\gcd(c_1,c_2)\).
Thus equal-order common target primes contribute only once.
For equal targets the charge is just their outside radical.

Write \(z_i^2D_i=a_i b_i\), and let
\(\beta_i=a_i-b_i+2z_i\sqrt{-D_i}\). Its norm is \(c_i^2\).
At odd \(p\mid c_i\), \(a_i b_i\) is a unit and
\(-D_i\) is a nonzero square modulo \(p\), since
\(-a_i b_i\equiv a_i^2\bmod p\). The two prime ideals
above \(p\) are unramified of degree one. The conjugates
of \(\beta_i\) sum to \(2(a_i-b_i)\), a unit there,
so exactly one of the two ideals is active, with
exponent \(2v_p(c_i)\). The other is the unit branch.
On this branch \(\beta_i\equiv4a_i\bmod
p^{v_p(c_i)}\), because its conjugate has at least that
valuation and \(a_i-b_i\equiv2a_i\bmod p^{v_p(c_i)}\).
Orient the active primes for each row uniformly over
the whole retained family. A random assignment of
the two choices at each relevant odd prime and row
keeps a fraction at least \(X^{-o_X(1)}\), using the
uniform bound on \(\omega(c_i)\). The number of relevant
primes is finite at each height.

Let \(C_*=\operatorname{lcm}(c_i)_i\) and remove from it all
primes dividing \(M\), writing the result \(W_*\).
Normalize the elements to \(\beta'_i=(W_*/w_i)\beta_i\).
These multipliers are integers and units at \(M\).
At the unit branches modulo \(M\), the residue vector
is that of \(4W_*(a_i/w_i)_i\). Its class modulo the
diagonal scalar is the stated element of \(S\).
There is no requirement on equality of target orders
on the primes of \(M\).

Put \(2M\) in the prescribed support, together with
the original \(P\). At \(p\) outside it list the distinct
positive target valuations \(0<k_1<\cdots<k_t\), with
\(k_0=0\). At the layer \(j\) take the tuple of ideals
with component \(\mathfrak p_{ip}^2\) if
\(v_p(c_i)\ge k_j\), and \((p)\) otherwise. Here
\(\mathfrak p_{ip}\) is the fixed active ideal, required
only in the active components. Raise this tuple to the
ordinary common exponent \(k_j-k_{j-1}\).
Its component norms are all \(p^2\), and at least
one component has proper support above \(p\).
The product over layers gives precisely the local
ideal tuple of \((\beta'_i)_i\). In coordinate \(i\)
it contains the active order \(2v_p(c_i)\) and the
rational factor of order \(\max_j v_p(c_j)-v_p(c_i)\).
The squarefree \(p\)-product in each slot records its
ideal tuple uniquely, once the nonempty active set and
the layer exponent are specified.

Fix the supported factors at height-subpower cost.
There are only a fixed number of layer occurrences on
each prime. The product of all layer primes with
multiplicities and exponents is at most \(C_*\le X^m\).
Thus the radical of the primes carrying any layer
exponent above \(J\) is at most \(X^{m/(J+1)}\).
Fix all the powers on these primes with a total price
less than \(\eta/8\) plus subpower by first choosing \(J\)
at each fixed \(L\). On the remaining slots use a finite
mesh with summed cap \(H\le\ell+\eta/4\). Each cap
\(n^\theta\) bounds the squarefree rational root and the
dictionary size, while the component norms are at most
\(n^{2\theta}\).

For each field let \(\mathfrak m_i\) be the modulus
made of its unit branches above \(M\), with the
exponents in \(M\). Use the product of ordinary ray
class groups modulo \(\mathfrak m_i\), divided by the
image of the diagonal rational unit residues. This is
a finite abelian group \(G_*\). Two ideal tuples
represent the same class exactly when their quotient
has generators with the same rational unit residue
in all the selected branches. The image of \(U\)
is the subgroup of principal ideal tuples modulo
that relation. Its kernel has size at most \(6^m\),
since each imaginary quadratic unit group has size
at most six. For each fixed positive \(e\), the ray
exact sequence and the local prime-power unit bound
give
\[
          |G_*[e]|\le X^{o_X(1)}
                         \prod_i |\mathrm{Cl}(K_i)[e]|.
\]
Passing to a quotient cannot increase \(e\)-torsion
size in a finite abelian group, since it equals the
order of the quotient by \(e\)-multiples. Minkowski's
bound and the ideal divisor count bound each class
number by \(O_\zeta(D_i^{1/2+\zeta})\) for fixed
\(\zeta>0\).

In a slot moment, choose the first product of \(h\)
roots. Among second products with the same \(e\)-multiple
in \(G_*\) there are at most \(|G_*[e]|\) possible
classes. For two second products in one of these
classes let \(\gamma_i\) generate their ideal
quotients, with identical rational residues at the
unit branches. All component norms of either
product are equal. If \(E\) is the common norm of
the denominator components, it clears all
denominators, and
\[
                |\sigma(E\gamma_i)|\le n^{2h\theta}.
\]
This is the imaginary quadratic norm identity for
the modulus of a complex embedding. If \(K_i,K_j\)
are distinct, work in their biquadratic compositum.
At every prime power \(p^a\parallel M\) there is a
degree-one prime above the chosen pair of branches.
There \(E(\gamma_i-\gamma_j)\) vanishes to order at
least \(a\). Hence a nonzero difference has integer
norm at least \(M\) in absolute value, whereas its
norm is at most \(16n^{8h\theta}\).
Choose \(h=\lceil T_0H/\theta\rceil\), with \(T_0\)
a fixed sufficiently large moment multiplier, and
take \(L\) so small that
\(8L(T_0+1)H<w\). Then the generators agree and
are rational, because the quadratic fields intersect
in \(\mathbb Q\). Compare each row with one in a
distinct field to obtain one common rational
generator in every component. Every template has
a proper component, which a nonzero rational
valuation would fill on both branches. Thus all
the valuations of that generator on the moving
primes vanish, and the second product is fixed.
Divisor-many ordered decompositions finish the
moment bound, with half the slot price and a
factor \(|G_*[e]|^{1/(2h)}\).

Apply Fourier inversion and Hölder with
\(\sum(2h)^{-1}\le1/(2T_0)\).
For the ordinary list, push its indicator to
\(G_*\) and use \(L^2\). For relative expansion,
use the RX1 majorant on \(U\), evaluated on
characters of \(G_*\). The kernel bound gives
the same \(L^1\) estimate up to a constant.
Choose \(\xi\) with \(b\xi<\eta/8\), then choose
\(T_0\) so that \(1/(2T_0)\le\xi\) and the
class number factors cost at most \(n^{\eta/8}\).
All these choices precede the small-budget threshold.
Each fixed final ideal has at most six generators.
A generator \(\beta'_i\) recovers its normalized
target as the positive square root of its norm,
and the difference of its normalized summands from
its real part. Recover the positive summands and
divide them by their gcd to obtain the unique
original primitive row. This completes the count.
[Ordinary ray class exact sequence, Stevenhagen, §3](https://websites.math.leidenuniv.nl/algebra/Stevenhagen-CFT.pdf)
· [Single absolute residue list in a ray group (AT1)](ABC_ABSOLUTE_TARGET_RAY_MOMENTS.md#algebraic-integers-with-a-ray-residue-dictionary)

### Unscaled absolute residues and balanced block prices

The count also applies with the absolute residue vector
\((a_i\bmod M)_i\), modulo the diagonal scalar, in place of
\((a_i/w_i\bmod M)_i\). Use \(\beta_i\) directly.
For each distinct positive order \(e\) at a prime, use the
template with component \(\mathfrak p_{ip}^2\) in the rows
with \(v_p(c_i)=e\) and the identity ideal in the other
rows. Its exponent is \(e\), so the same distinct-order
radical charge applies. The moment comparison now clears
denominators by the product of the denominator norms, with
\(m\) fixed. If the slot cap is \(n^\theta\), each conjugate
after clearance is bounded by \(n^{(2m+1)h\theta}\).
Taking \(L\) smaller in terms of \(m,T_0,H,w\) preserves
the norm contradiction. No assumption of equal orders
at the primes of \(M\) is needed. In particular a genuine
absolute residue cell supplies the singleton class in the
diagonal quotient. Neither a catalogue of its absolute
residue value nor an equality of complete targets is
required for the joint bound.

**RX4.** In the setting of the preceding unscaled count,
let \(\ell_1,\ell_2\) cap the two original outside target
radicals. Suppose eligible ordered pairs have distinct
parity fields and an actual field-pair catalogue of size
at most \(n^{f+o_X(1)}\). Suppose the outside radical
mass of common target primes of equal positive order
is at least \(\chi\). Let \(\zeta\) be a lower bound on
the radical mass of the target primes that occur in the
other row's parity field kernel, outside the prescribed
pool and modulus. Such a prime occurs in just one target
by original primitivity. The two overlap sets are disjoint.
Include both prescribed field supports in the pool.
They remove the mass \(\zeta\), and common equal-order
templates remove \(\chi\) from the sum of the original
charges. In a common residue class the total pair
exponent is at most
\[
             f+\frac{\ell_1+\ell_2-\chi-\zeta}{2}+\eta .
\]
If any pool records had to be fixed with catalogue
price \(b_1\), add \(b_1\). The union of the two
squarefree field kernels outside the initial pool
furnishes a sufficient actual field-record price:
fix the union at its radical cost and allocate each
prime to either or both kernels at subpower cost.

For example, in a retained AP2 block put \(c_0=1/6\),
let \(\beta\) be the normalized radical mass lower bound
for its target divisor, and let \(s\) be its packet
price as in AP2. The odd part of the packet
supplies a fixed positive numerical lock, because
its order at \(2\) is bounded at each fixed \(L\)
by AP1. Use the known full packet support as part of \(P\).
A block of size \(n^{2(c_0-s)+o(1)}\) with an
\(n^{-\rho+o(1)}\) fraction of the indicated ordered
pairs violates the bound if
\[
               b_1+f <
                  2c_0-4s+\beta+(\chi+\zeta)/2-\rho
\]
by a fixed margin. A fixed positive block tolerance
only changes the formula by the stated AP2 losses
and must be chosen inside that margin. This follows
by squaring the block size and using
\(\ell_i\le2c_0-\beta+o(1)\).
At the independent field price \(2(c_0-s)\), with
\(b_1=0\), a strict saving requires
\(\chi+\zeta>4s-2\beta+2\rho\). Choosing only
equal-target pairs in an evenly occupied block
pays \(\rho=c_0-s\). If both raw target caps are
saturated and there is no supplementary record
compression, the full target overlap available
there gives no strict saving, since \(s\ge\beta/2\).
The useful missing input is a joint field catalogue
or overlap exceeding the stated extraction cost.
The estimate compares genuine pairs with their
original three caps throughout.
[Composite residue packets (AP2)](ABC_AGGREGATE_DEEP_PACKET_INVERSE.md#composite-packet-inverse)

### The cost of extracting one-directional overlaps

**RX5.** In the AP2 block above, write \(d=c_0-s\).
Separate the cross-field mass \(\zeta\) into a mass
\(\zeta_{12}\) of primes from the first row's field
kernel on the second target, and \(\zeta_{21}\) in
the other direction. All mass parameters are in a
fixed fine mesh. If the pair sector has exponent
\(4d-\rho+o(1)\), then
\[
         2\rho\ge\chi+\max(\zeta_{12},\zeta_{21})
                              -2s+\beta-o(1).                \tag{RX5}
\]
In fact with the first row fixed, include the
support of \(c_1 D_1\) in the second target's
prescribed pool. This is legitimate even when
the original field catalogue is large, because
the actual first row has already been fixed.
The remaining second target radical is at most
\(n^{2c_0-\beta-\chi-\zeta_{12}+o(1)}\).
The two overlap sets on that target are disjoint
by primitivity of the first row. Fix one complete
summand of the second row. FL1, with this pool
on its target, bounds the completion count by
half the remaining target exponent, with an
arbitrarily small fixed allowance. There are
\(n^{d+o(1)}\) different values of the chosen
summand in the actual block. Multiplying by
the first-row count \(n^{2d+o(1)}\) proves the
one-directional bound for \(\rho\). Exchanging
rows and taking the larger restriction proves
RX5. If a genuine near-full block has actual
outside target radical cap below \(2d\) by a
fixed margin, the same fixed-summand count
already contradicts its size.

Consequently at the independent field-pair price
\(2d\), target overlap alone physically cannot
finance the strict saving in RX4, even across
different targets. For those inequalities to
deliver a strict saving without a field catalogue
discount or a smaller raw cap, one needs
\[
                    \min(\zeta_{12},\zeta_{21})>2s-\beta .
\]
Indeed
\(\chi+\zeta-2\rho\le
\min(\zeta_{12},\zeta_{21})+2s-\beta+o(1)\)
by RX5, while RX4 needs
\(\chi+\zeta-2\rho>4s-2\beta\).
This compares exactly the charges furnished by
the two estimates. It identifies the need for a
joint field discount or target-to-operand overlaps
in both directions at that catalogue price.
For a fixed field pair there is also an anchored
comparison. If \(r\) caps each original outside
target mass in the block, the first row has at
most \(n^{d+o(1)}\) choices at its fixed field:
there are \(n^{d+o(1)}\) complete targets, and for
each fixed target its quadratic norm has only
divisor-many ideals and at most six generators
per ideal. Given the first row, the second row
has exponent \((r-\chi-\zeta_{12})/2\) by AT4
with its fixed field, fixed absolute residue and
the same enlarged pool. Splitting its prescribed
squarefree parity kernel into the two coprime
operand kernels costs only a subpower. Interchange
the roles to get the summed anchored exponent
\[
           f+d+\frac{r-\chi-\max(\zeta_{12},\zeta_{21})}{2}.
\]
Thus the joint ray exponent improves this anchored
exponent exactly when \(\min(\zeta_{12},\zeta_{21})>
r-2d\), up to arbitrary fixed allowances and with
the same actual field catalogue. The overlap must
use both original target-to-operand directions to
give an advantage over this comparison at
\(r=2c_0-\beta\).

If the only available field price comes from
enumerating the actual kernels by their radical,
there is a stronger anchored bound. Let \(f_0\)
cap the normalized field union outside the initial
pool and both actual targets. Stabilize the two
cross masses in a fixed fine mesh. The raw field
price is then \(f_0+\zeta_{12}+\zeta_{21}\).
Fix the first complete target before fixing the
fields, at cost \(d\). The mass \(\zeta_{21}\)
now lies on a prescribed support and costs a
subpower. The remaining field union costs
\(f_0+\zeta_{12}\), then the first row is determined
up to subpower by its target and field. Count the
second row as above and choose the better orientation.
The resulting bound has exponent
\[
       f_0+d+\frac{r-\chi+\min(\zeta_{12},\zeta_{21})}{2}.
\]
This is at most the raw-price joint ray exponent
\(f_0+r-\chi/2+(\zeta_{12}+\zeta_{21})/2\) in a
near-full block, where \(r\ge2d-o(1)\). Thus a
genuine joint catalogue discount is needed to improve
this anchored comparison when raw support enumeration
is the only field input.
[Fixed-letter counts with prescribed pools (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)
· [Single absolute input residue (AT4)](ABC_ABSOLUTE_TARGET_RAY_MOMENTS.md#divisors-of-binary-norms-with-an-absolute-input-residue)

### Finite ratio templates and the same-field exception

A fixed finite catalogue of primitive valuation vectors
can be charged once per prime. At a prime with
\(v_p(c_i)=t r_i\), where \((r_i)_i\) is in the
catalogue and \(r_*=\max_i r_i\), replace the layer
templates by the one template with component
\(\mathfrak p_{ip}^{2r_i}(p)^{r_*-r_i}\) at
common exponent \(t\). An absent active factor
means the identity. All component norms are
\(p^{2r_*}\), and a maximal coordinate is proper.
This merges primes in the catalogue into one unit
of charge. Use ordinary layers at all other primes.
The moment argument goes through with its numerical
threshold adjusted by the fixed bound on \(r_*\).
This catalogue is fixed before choosing the small
budget threshold. A list growing to accommodate
all valuations is not covered by that argument.

In a common quadratic field, the common-residue
diagonal can allow an equal nonrational quotient.
For example take \((\gamma_1,\gamma_2)=((2+i)^2,
(2+i)^2)\) in \(\mathbb Q(i)\) on the same unit
branch above \(13\), for arbitrary \(M=13^k\).
The ideal tuple is supported on a proper branch
above \(5\) in both components. This is an exact
counterexample to the rational-generator step alone.

For two rows in the same field, one can still
obtain the bound if the full parts of \(M\) where
the prescribed unit branches align and where they
are conjugate both have a fixed positive numerical
weight. The two comparisons then force
\(\gamma_1=\gamma_2=\overline{\gamma_2}\).
There is also a uniform alternative without this
split. At least one alignment part has numerical
weight at least \(w/2\). Conjugate one coordinate
globally if necessary, and the comparison on that
part forces \(\gamma_1=\gamma_2\).
Fix at full radical price the common layers on which
both active branches agree in this orientation.
Their radical counts each such prime just once
since in the ordinary two-row layer templates
each prime has at most one common layer. In a
finite mesh for its actual normalized radical mass
\(\kappa'\le\kappa\), record those support primes
and their aligned common layers at price
\(\kappa'\) up to the mesh allowance. The other
layers have charge at most \(\ell-\kappa'\).
In the other slots equality of the two ideal quotients
makes every valuation zero. For an exclusive
layer one coordinate has a rational ideal while
the other is proper. For a common layer of
opposite branches their ideal supports are
different. The bound is therefore \(n^{\ell/2+
\kappa/2+\eta}X^\varepsilon\) times the same list
factor, where \(n^\kappa\) bounds the radical of
these fully aligned common-layer primes. This
alternative uses the ordinary layer price. It
does not require distinct rows, and with disjoint
outside target radicals the extra cost is zero.


## A two-cap example with maximal normalized word growth

**RX6.** Fix \(0<f<1/6\), \(B_*>0\), and a finite list
of signed integer words, each with fixed nonzero slot
weights. For arbitrarily large fixed even \(k\), put
\(L=3/k\) and \(n=X^L\). For all sufficiently large
inner heights there is a matching of
\(n^{1/3-f-o_X(1)}\) distinct positive primitive
rows \(A+B=C\) with the following properties.
All three letters have fixed positive numerical
height exponent. Both \(\operatorname{rad}(A)\)
and \(\operatorname{rad}(C)\) are at most \(n^{1/3}\).
All targets share one full prime power divisor
\(D\) with \(\operatorname{rad}(D)\le n^f\) and
\(D\ge X^{2f}\). Every row has
\(\operatorname{rad}(B)>n^{B_*}\).
For the normalized residues \(\theta=A/(C/D)\bmod D\),
each prescribed word of length \(j\) has image size
\(n^{j(1/3-f)-o_X(1)}\) on ordered tuples of
distinct retained rows. This supplies maximal word
growth with the actual target cap as well as one
operand cap.

Choose an odd prime \(\pi\) in \([n^f/2,n^f]\)
and \(D=\pi^k\). For primes \(p,q\) in the fixed
relative intervals giving
\[
          A=p^k\in[X/8,X/4],\qquad
          C=(\pi q)^k\in(X/2,X],
\]
put \(B=C-A\). The prime number theorem in intervals
of fixed relative length gives \(n^{1/3-o_X(1)}\)
choices of \(p\) and \(n^{1/3-f-o_X(1)}\) choices
of \(q\), uniformly in the chosen \(\pi\).
Their scales are distinct from each other and
from the scale of \(\pi\), because \(f<1/6\).
This makes the rows primitive, and \(D\) is the full
\(\pi\)-part of \(C\). In particular \(B\) is also
bounded above and below by constant multiples of
\(X\).

To discard small radical on \(B\), use the
\(\tau(k)\) homogeneous cyclotomic factors of
\((\pi q)^k-p^k\). Off the fixed primes of \(k\)
their supports are disjoint, by separability and
the coprimality of \(\pi q,p\). If
\(\operatorname{rad}(B)\le n^{B_*}\), at least one
of these factors has outside radical at most
\(n^{B_*/\tau(k)}\). For fixed \(p\), each such
value costs at most \(n^{B_*/\tau(k)+o_X(1)}\)
choices by supported-power enumeration, and
the nonconstant polynomial in \(q\) has at most
its degree many roots for a specified exact value.
Take \(k\) through multiples of a sufficiently
composite fixed even integer so that
\(B_*/\tau(k)<1/3-f\). This discards a vanishing
fraction of edges uniformly at each \(p\).
An edge conflicts with at most
\(O(|\{p\}|+|\{q\}|)\) other edges by a repeated
\(A\), \(C\), or \(B\); for the last claim each
\(q\) permits at most one \(p\). Greedy selection
therefore keeps \(n^{1/3-f-o_X(1)}\) edges with
all three complete dictionaries injective.

On the matching \(\theta=(p/q)^k\bmod D\).
For any one fixed weighted word on distinct rows,
the reduced rational word in the fractions \(p/q\)
has height at most \(X^{\sum_i|e_i|/k}\).
Each prime belongs to just one matched index,
and the exact fraction determines the chosen
distinct indices and their weights up to finitely
many permutations. A specified word output modulo
\(D\) has at most \(2k\) possible \(k\)-th root unit
classes. Choose \(k\) so large that the square of
the rational height bound is smaller than \(D\).
Rational cross multiplication makes each root
class contain at most one of the exact fractions.
Thus each modular output has bounded
multiplicity and the stated image count follows.
Choosing one large enough \(k\) works for the
fixed finite list of words. The missing operand
cap is therefore essential to any argument
excluding this type of modular expansion.
[Prime number theorem, Soprounov, Theorem 2.1 at modulus one](https://academic.csuohio.edu/soprunov-ivan/wp-content/uploads/sites/93/2023/02/primes.pdf)
· [Fixed-support counting (EB3)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)



### Diffuse profiles and adaptive full target locks

**RX8.** The two-cap model also admits wholly diffuse scaled
valuation profiles and arbitrary adaptive choices of full
target divisors. Fix \(0<f<1/6\), a missing-operand lower
radical threshold \(B_*>0\), a positive numerical lower bound
\(\nu\) for the tested locks, and a finite list of fixed
nonzero signed words. Take arbitrarily large even \(k\)
with sufficiently many divisors and put \(L=k^{-2}\).
There are actual primitive matchings of size
\(n^{1/3-f-o_X(1)}\) with both the \(A\) and \(C\) caps
as in RX6 and with every \(\operatorname{rad}(B)>n^{B_*}\).
All numeric heights are comparable to \(X\).
In the outer limit, each of the two capped valuation
measures tends to \((1/3)\operatorname{Unif}[2,4]\).
On every fibre of every fixed actual full divisor
\(M\mid C,\ M\ge X^\nu\), each tested signed word in
\(A/(C/M)\bmod M\) has bounded multiplicity on
ordered tuples of distinct rows, uniformly over the
locks at any fixed \(k\). This remains true for
every subset of each fibre.

Use slots \(j=2k+1,\ldots,4k-1\), put \(J=2k-1\)
and use exact exponents \(e_j=kj\).
Their scaled average is three. Use disjoint prime
alphabets in each slot at the three respective
scales \(n^{1/(3J)}, n^{f/J},
n^{(1/3-f)/J}\). The first and third
are the varying numerator and residual target
primes \(q_j,r_j\). Choose one fixed prime
\(\pi_j\) per slot at the middle scale and
put
\[
 s=\prod_j q_j^j,\quad r=\prod_j r_j^j,\quad
 D=(\prod_j\pi_j^j)^k,\qquad A=s^k,\ C=Dr^k .
\]
These disjoint prime pools give primitive rows.
Write \(E_0=\sum_j e_j\). On the first scale
restrict the prime size relative to its scale
to \([8^{-1/E_0},4^{-1/E_0}]\). On the
middle scale use \([1/2,e^{\tau/E_0}/2]\)
and on the third \([2e^{-\zeta_2/E_0},
2e^{-\zeta_1/E_0}]\) for fixed
\(0<\tau<\zeta_1<\zeta_2<\log2\).
The prime number theorem gives disjoint
pools of the required exponents. The intervals
ensure \(A\in[X/8,X/4]\), \(C\in(X/2,X]\)
and \(B\ge X/4\). The \(A\) radical lies
below its raw cap, the \(D\) radical is
at most \(n^f\), and at the combined target
scale its prime multipliers multiply to
at most \(e^{(\tau-\zeta_1)J/E_0}<1\).
Thus \(C\) meets its cap. Both measures have
mass \(1/(3J)+o_X(1)\) in each scaled slot.
The displayed weak limit follows.

Let \(N_A,N_C\) count the two full dictionaries
of exponents \(1/3\) and \(1/3-f\).
Each Cartesian choice gives one edge, and
RX6's cyclotomic argument removes at most
\(N_A n^{B_*/\tau(k)+o_X(1)}\) bad edges,
where \(\tau(k)\) here denotes the divisor
count of \(k\). Choose \(\tau(k)\) large
enough by taking multiples of a fixed
integer. Independently sample edges with
probability \(1/(N_A\log X)\). The expected
size is \(N_C/\log X\). Repeated complete
\(A\), \(B\), or \(C\) labels and removed
edges number \(o(N_C/\log X)\) in
expectation, since a label occurs at most
\(\max(N_A,N_C)\) times.

For any nontrivial exact multiplicative
relation in the fractions \(s/r\) on \(t\)
distinct sampled labels, every prime in
every individual slot occurs on at least
two labels. This follows from unique
factorization, since the weights are
nonzero and the slot alphabets are disjoint.
There are at most
\(O_k((N_A N_C)^{\lfloor t/2\rfloor})\)
candidate label lists, by first choosing
their boundedly many equality patterns.
Each list of distinct edges survives with
the \(t\)-th power of the sampling
probability. The expected number of
such relations tends to zero, uniformly
for the finite collection obtained by
combining two tested words and merging
equal labels. Thus a sample of the
stated order can be chosen with none of
these relations and with only a vanishing
fraction of repeated or bad edges.
Delete those edges. The exact tested
words then have only the bounded
multiplicity due to permutations of labels.

On a fibre of an actual full divisor \(M\)
the \(k\)-th root of \(A/(C/M)\) differs
from \(s/r\) by one fixed rational scalar.
Here every prime power in \(C\) is a
multiple of \(k\), so the root is rational
and is a unit at \(M\). It has height
at most \(X^{1/k}\). For each tested word,
the rational root output therefore has
height at most \(X^{\sum_i|e_i'|/k}\)
when its weights are \(e_i'\).
Take \(k\) large in terms of the prescribed
weights and \(\nu\) so its squared height
is below \(M\). Each fixed modular \(k\)-th
power has at most \(k^{2J}\) unit root
classes, since \(M\) has at most \(2J\)
primes and all exceed \(k\) at large \(X\).
One root class determines the reduced
rational root output by cross multiplication.
The bound applies to every such fibre
and all its subsets with the same
constants. Thus diffuseness, actual
target eligibility and adaptive full
locks are compatible with maximal
normalized word spread under the
two stipulated original caps.
[Prime number theorem at modulus one](https://academic.csuohio.edu/soprunov-ivan/wp-content/uploads/sites/93/2023/02/primes.pdf)
