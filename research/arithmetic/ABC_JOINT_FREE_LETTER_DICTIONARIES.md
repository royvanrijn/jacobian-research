# Joint dictionaries across the free letters and moving quadratic kernels

A full divisor of one letter gives a multiplicative congruence on the
signed ratio of the two free letters. Fourier moments can use a
dictionary of coordinates from **both** free letters at once. At every
balanced limiting profile with a cheap positive-weight modulus on the
third letter, restricted joint dictionaries of this kind have full
entropy to first order. This also lets a first-moment theorem for
quadratic class groups replace a higher-moment hypothesis in a
moving-kernel sector of the mixed signature \((6,2,3^a)\).

## A signed joint dictionary estimate

Use the primitive anatomic boxes of (GD3). The budget \(L\) is
fixed before the height limit. For the free letters \(U,V\), write
\(\alpha_i\) for the raw root exponents in base \(X\), with
\(\sum_i\alpha_i=O(L)\), and \(e_i\) for their exact original degrees.
Suppose the remaining letter \(T\) has a prescribed block of full
prime powers with combined numerical height at least \(X^w\),
for a fixed \(w>0\), and all its prime exponents at least
\(\kappa/L\), for fixed \(\kappa>0\). A record of the block values
costs at most \(X^{\mu L+o_X(1)}\) possibilities. Fix the
exceptional coefficients at the ordinary truncation loss.

Partition the free root positions into groups \(G\) that may meet
both \(U\) and \(V\). Each group is allowed boundedly many positions,
or its scaled degrees \(Le_i\) lie in a bounded number of windows of
diameter tending to zero with \(L\). The size and window-count caps
are fixed in the small-budget limit. Let the joint dictionary on
\(G\) have cardinality at most \(X^{W_G+o_X(1)}\), with
\(W_G\le\sum_{i\in G}\alpha_i\). One extra group, with no bound on
the number or locations of its degrees, is allowed if
\(\sum_{i\in G}e_i\alpha_i\le w-\xi\) for some fixed \(\xi>0\).
For every fixed \(\eta>0\) and all sufficiently small fixed budgets,
the count is at most
\[
            X^{\mu L+\frac12\sum_G W_G+\eta L+o_X(1)}.
 \tag{JD1}
\]
One can allow several exceptional groups by choosing positive
integer moment orders \(h_G\) for them with
\(h_G\sum_{i\in G}e_i\alpha_i \le w-\xi\), provided the moment
weights described below still sum to at most one.

[Grouped moments within a letter](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md#diagonal-moments-for-a-joint-dictionary)
· [Signed monomial congruences (FC1–FC2)](ABC_OMITTED_COORDINATE_GRAPH.md#signed-monomial-congruences-over-full-moduli)
· [The wider-cofactor bound (GW)](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-bounded-number-of-small-roots)

To prove (JD1), color each prime by a single original position,
including the letter in the position label. Independent random
coloring retains \(X^{-o_X(1)}\) of all rows with each prime in
its assigned position, since the number of positions is fixed at
each budget and \(\omega(abc)=O(\log X/\log\log X)\). Restrict each
dictionary to the coloring. At a fixed modulus \(M\) from the block
on \(T\), every retained row satisfies
\[
       \prod_G r_G\equiv K\pmod M,\qquad
       r_G=\prod_{i\in G}x_i^{\sigma_i e_i},
       \qquad \sigma_i=1\text{ on }U,\quad \sigma_i=-1\text{ on }V .
\]
The unit \(K\) incorporates the fixed exceptional coefficients and
the signs of the additive equation. Within each admissible group
dictionary discard factors that are not units modulo \(M\). Let
\(S_G(\chi)=\sum_{\mathbf x_G}\chi(r_G)\) for a group dictionary on
the character group of \((\mathbf Z/M\mathbf Z)^\times\). Count all
joint dictionary choices with the displayed residue by
\[
     \mathbb E_\chi \overline{\chi(K)}\prod_G S_G(\chi).
\]
The exact choices determine \(U,V\), so each determines at most one
original triple.

Write \(\alpha_G=\sum_{i\in G}\alpha_i\) and
\(L_{UV}=\sum_G\alpha_G\). Use
\(h_G=\lceil L_{UV}/\alpha_G\rceil\) on ordinary groups of positive
cap. Zero-cap groups have subpower dictionaries and can be fixed.
Upper root caps and interval refinements can absorb arbitrarily small
total rounding loss at each fixed budget before choosing the
moment orders. For each ordinary group, expand
\(\mathbb E_\chi|S_G(\chi)|^{2h_G}\) and fix the first \(h_G\)
tuples. Collapse the second \(h_G\) tuples coordinatewise. The
unsigned product in each group has height at most
\(X^{2L_{UV}+o_X(1)}\). A collapsed value has divisor-many ordered
preimages, at fixed truncation and moment order. Primes from distinct
original positions remain disjoint by the coloring. For a bounded
group, (FC2) counts the signed congruence on the collapsed roots at
\(X^{\eta L+o_X(1)}\) cost. Its amplification order can be chosen
large before \(L\) becomes small, uniformly in the group.

For a bounded number of shrinking windows, compress the collapsed
roots to one root per window and per sign. Use a common degree
representative in each positive window and treat the exponent
differences, together with any window tending to zero, as a signed
wider cofactor of height \(X^{o_{L\to0}(1)+o_X(1)}\). This follows
because \(h_G\alpha_G=O(L)\). The remaining positive representative
degrees must be bounded below after multiplying by \(L\), as required
in (GW). To arrange this for varying centers, pass to a subsequence
in the budget limit on which each of the boundedly many scaled
centers in the group tends to zero, to a positive scale, or to
infinity. A window tending to infinity has representative degree
at least a fixed multiple of \(1/L\) as well. The same argument
applies uniformly to any sequence of groups witnessing a failure of
the moment estimate. The signed wider-cofactor bound gives
\(X^{\eta L+o_X(1)}\) possibilities, using its amplification at fixed
order and the narrow root height cap. The exact signed monomial
recovers all original products from their prime colors and exact
degrees. This also controls reconstruction where the wider
cofactor and grouped roots use the same support. Thus an ordinary
moment has bound \(X^{h_G W_G+\eta L+o_X(1)}\).

For an exceptional group, use \(h_G=1\), or its prescribed larger
order. If two \(h_G\)-fold signed group products agree modulo \(M\),
the cross-difference of their reduced numerators and denominators
has absolute value at most
\(2X^{h_G\sum_{i\in G}e_i\alpha_i+o_X(1)}<M\).
It is zero by divisibility. The coloring recovers the collapsed
products as before, giving the same moment bound with just
\(o_X(1)\) loss. Finally \(\sum_{\rm ordinary}1/h_G\le1\), and
one exceptional group with \(h_G=1\) gives
\(\sum_G1/(2h_G)\le1\). More generally impose this last condition
on the moment orders. Apply Hölder with exponents \(2h_G\),
using a constant function for any unused weight. Each group
contributes \(W_G/2\), and the moment errors are averaged with
total weight at most one. Sum over the modulus records.

### Consequences at the balanced endpoint

In the iterated limit suppose the count exponent divided by \(L\)
tends to \(1/3\), the raw radical mass on each letter tends to
\(1/3\), and \(T\) has positive-weight full moduli at arbitrarily
small normalized cost. Then for any allowed partition of
the \(U,V\) coordinates into ordinary groups,
\[
       \sum_G H(\mathbf x_G)
             =\left(\frac23+o(1)\right)L\log X .
 \tag{JD2}
\]
The statement also holds with one exceptional group whose total
numerical weight tends to zero in the small-budget limit. All
entropies are for the uniform law on the primitive rows in the box.
Indeed, bin the realized group surprises in a mesh with arbitrarily
small total normalized error at each fixed truncation, capping by
the raw group sizes when necessary. The number of profiles is
independent of \(X\). Any profile with a fixed cost deficit has a
power-smaller count by (JD1), after taking \(\mu\) and \(\eta\)
smaller than that deficit. In the exceptional case fix the cheap
positive-height block before going to budgets where the group
height tends to zero. Averaging the capped costs gives the lower
bound on the uncapped entropies, and the raw caps give the upper
bound.

For example one can match original positions across \(U,V\) in
any way and sum the pair mutual informations. Their sum is
\(o(L\log X)\), by applying (JD2) to the matched pairs and singleton
remainders and bounding each marginal by its raw cap. The same
proof works with arbitrary bounded-window matched blocks and on
any near-saturating subfamily. If a fixed record with surprise
\(s\) is first observed, the same binning gives the conditional
version with \(H(\text{record})\) added to the right of the
entropy upper bound from (JD1). The record need only give group
dictionaries in each of its fibers. No assertion of full entropy
on the unrestricted joint pair follows. In particular these
bounds apply beyond the minimum nonatomic mass profile.
[Entropy at minimum nonatomic mass (EN3)](ABC_DIFFUSE_FACTOR_GRAPH.md#consequences-of-the-saturated-count)
· [Cheap full moduli and zero nonatomic alternatives](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees)

## An unconditional sector with fully moving square kernels

Fix \(a\ge2\), \(0<F<1/3\), a positive integer \(m\), and an
orientation of \((6,2,3^a)\). Use the notation \(U,V,T\) in
that order, with all three complete values positive and the
appropriate sign for addition. Call \(S\) the squarefree product
of primes whose complete valuation on its own letter is not
divisible by that letter's order \(6,2,3^a\). On \(U,V\)
call the full parity-defect product
\[
           B_{\rm odd}
              =\prod_{\substack{p^e\parallel UV\\e\text{ odd}}}p^e .
\]
The canonical squarefree kernels \(A,B\) of \(U,V\) are
coprime, and \(AB=\operatorname{rad}(B_{\rm odd})\mid S\).
There are constants \(\delta,\gamma,\tau,\lambda_0>0\), depending
only on \(a,F,m\), with the following property. For every fixed
\(0<L<\lambda_0\), count primitive triples \(a_1+b_1=c\le X\)
with \(\operatorname{rad}(a_1b_1c)\le c^L\) and \(S\le c^{FL}\).
Impose at least one of these two conditions:
* \(B_{\rm odd}\le c^\gamma\).
* All but \(\tau L\log c\) of the logarithmic radical of
  \(B_{\rm odd}\) occurs at odd degrees with \(Le\) in a union of
  at most \(m\) real intervals of length at most \(\gamma\).
  The intervals may vary with the row.

Then the count is
\[
          O_{a,F,m,L,\varepsilon}
                  (X^{(1/3-\delta)L+\varepsilon}).
 \tag{JD3}
\]
The parity kernels and all intersection patterns with \(S\) may vary.
The same estimate counts the primitive quotients of least repairs.
For the balanced reverse transfer, if \(H_R(c)\le Y\) and
\(\operatorname{rad}(c)\ge c^{(1/3-\xi)L}\), the exponent becomes
\((1/3-\delta)/[R(1/3-\xi)]+\varepsilon\), strictly less than
\(1/R+\varepsilon\) at sufficiently small fixed \(\xi\).

### Torsion splitting and first-moment compression

Put \(n=X^L\). Work first in a dyadic height shell, so restrictions
using \(c\) imply the indicated \(X\)-caps. The cyclic isogeny
descent for \((6,2,3^a)\) gives a constant \(C_a\) with at most
\[
          \exp O_a(1+\omega(S))\,
             \bigl(h_3(\kappa)h_3(-3\kappa)\bigr)^{C_a}
 \tag{JD4}
\]
solutions per \(A,B,S\), even with no height bound. Here
\(\kappa=\pm AB\) according to the signed orientation, and
\(h_3(d)=|\operatorname{Cl}(\mathbf Q(\sqrt d))[3]|\),
with the trivial field assigned value one. This is (RR5) with
the cyclic isogeny and the quadratic class decomposition. For
each \(S\) the allocated kernels have at most \(3^{\omega(S)}\)
possibilities. Hence the rows on which the torsion product in
(JD4) is at most \(n^\theta\) have count
\(n^{F+C_a\theta+o_X(1)}\). Fix \(\theta>0\) so small that
\(F+C_a\theta<1/3\), and take \(\tau<\theta/16\).
[Isogeny descent and reduction to quadratic torsion](ABC_ROOT_RATIO_MIXED_COVERS.md#cyclic-isogeny-descent-for-2p2pa)

For the other rows, every discriminant with \(AB\le n^f\)
has at least one of the two class numbers \(h_3(\kappa)\),
\(h_3(-3\kappa)\) greater than \(n^{\theta/2}\). Both quadratic
discriminants have magnitude \(O(n^f)\), and the parameter \(AB\)
is recovered with bounded multiplicity from either discriminant.
The Davenport–Heilbronn first moment and Markov's inequality give
at most
\[
                        O(n^{f-\theta/2})
 \tag{JD5}
\]
possible squareclasses, including no such class once this
bound falls below one. Kernel allocations and the full original
valuations on their support have subpower multiplicity, uniformly
at polynomial height. Thus (JD5) bounds the full joint dictionary
of the odd exact-degree roots, up to the inner subpower loss.
[Quadratic \(3\)-torsion first moment, Davenport–Heilbronn, Theorem 3](https://www.math.utoronto.ca/~ila/DH2.pdf)

### Excluding high torsion at all relevant balanced profiles

Suppose no uniform saving (JD3) holds. With \(\tau\) fixed as
above, take candidate values of \(\delta,\gamma,\lambda_0\)
tending to zero and choose violating fixed budgets \(L\to0\),
with the height limit taken at each such budget. Split into the
two conditions and the low and high torsion sectors. The latter
is required by (JD4). Dyadic anatomic refinement, the limiting
one-third bound and (GD3) give a sequence of balanced boxes
where the normalized count tends to \(1/3\) and each letter's
radical mass tends to \(1/3\). In particular the sum of actual
root masses divided by \(L\) tends to one, since otherwise one
could use the limiting bound at a uniformly smaller budget.
All prime powers on the support of any fixed exceptional factor
can be enumerated at subpower cost. In the standard fixed
truncation their radical loss can be taken arbitrarily small
relative to \(\theta L\). Equivalently, fix their full records
and remove them from the group to which (JD5) will be applied.
An exceptional tail on the modulus letter of positive numerical
weight itself supplies a full modulus, with the radical cost
bounded by the truncation loss.
[Limiting one-third count](ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13)
· [Balance from the two-letter bound](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md#group-dictionary-count)

Refine so the combined odd kernel raw cap is \(n^{f+o_X(1)}\),
taking the raw mass as \(f\). In the window case, the window
classifications of exact degrees have only finitely many possible
records at fixed truncation, so choose a largest class. The group
of odd root positions in those windows has raw mass at least
\(f-\tau-o(1)\), in normalized units. Its joint dictionary is
a projection of (JD5) and has cap at most \(f-\theta/2+o(1)\).
In the small-height case use the whole odd group. Its numerical
weight tends to zero; it has raw mass \(f+o(1)\) and the same
dictionary cap. All other positions can be singletons.

If \(T\) has arbitrarily cheap positive-weight full moduli along
the limiting sequence, choose one at cost \(\mu<\theta/32\).
In the small-height case the odd group qualifies as an exceptional
group in (JD1) because its numerical weight tends to zero.
The group in the window case is ordinary. In either case the
pair has a total dictionary deficit of at least
\(\theta/2-\tau-o(1)\). Equation (JD1) now gives a fixed
gap below \(1/3\), a contradiction.

If \(T\) lacks such moduli, the zero-height and finite-scale
alternatives at balance make the two other limiting degree
measures wholly nonatomic. A bounded number of shrinking windows
on \(U,V\) then carries mass \(o(1)\). Indeed compactness gives
a subsequential limit for every bounded window center and no
atom can occur there. The bounded numerical first moments give
vanishing radical mass for centers tending to infinity. In
the windowed sector this proves \(f\le\tau+o(1)\). In the
small-height sector all odd mass at scaled degrees above
any fixed \(t>0\) is bounded by the tending-to-zero odd numerical
weight divided by \(t\). The mass below \(t\) tends to zero
as \(t\downarrow0\), by the absence of atoms at zero. Here
\(f=o(1)\). Both conclusions contradict the existence of any
high-torsion class in (JD5) for late budgets. This proves the
existence of the uniform constants in (JD3). Summing over
height shells or applying the unique least repair preserves
the stated exponent allowances.
[Finite-scale and zero-height alternatives](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees)
· [Least-target repair](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)

### A necessary odd-degree mass at saturation

Fix \(\theta>0\) with \(F+C_a\theta<1/3\). At a balanced limiting
endpoint in this signature-defect sector, form the two submeasures
on scaled degrees retaining only the exact odd valuations on
\(U,V\). Their combined nonatomic mass must be at least
\(\theta/2\). If the combined nonatomic mass were smaller, choose
finitely many atoms, including zero when it has positive mass,
leaving outside mass strictly below \(\theta/2\). Slowly shrinking
continuity windows around those atoms capture the chosen mass in
the iterated limit. The raw odd-root group then saves the difference
between \(\theta/2\) and the remaining mass in (JD5). Choose the lock
cost smaller than half that saving and use (JD1). If no such lock is
available, the wholly nonatomic measures on \(U,V\) give
\(f<\theta/2\), which rules out the high-torsion sector directly.
The cutoff at large scaled degrees costs arbitrarily small mass
by the uniform first-moment bound.

### The unrestricted kernel group in this mixed signature

For fully diffuse odd-kernel degrees of positive numerical weight,
a sparse list of total parity kernels need not give a sparse
bounded-window projection. The unresolved moment can be stated
directly. Use the same high-torsion threshold \(n^\theta\) and the
raw odd mass \(f\) in a refined box. For a fixed cheap full modulus
\(M\) on \(T\), let \(\mathcal D_{\rm odd}(M)\) be the joint odd
dictionary from (JD5), restricted to units modulo \(M\). Write
\(r_{\rm odd}\) for its full signed prime-power ratio. A sufficient
bound on this dictionary, for some fixed \(\sigma>0\), is
\[
 \sum_{b\in(\mathbf Z/M\mathbf Z)^\times}
       \big|\{\mathbf x\in\mathcal D_{\rm odd}(M):
                   r_{\rm odd}(\mathbf x)\equiv b\bmod M\}\big|^2
       \ \ll\ n^{f-\sigma+o(1)} .
 \tag{JD6}
\]
Here uniformity in \(M\) at arbitrarily small positive locking cost
suffices, with a fixed positive limit of the numerical weight at
each such cost. By taking \(h_{\rm odd}=1\) and keeping the ordinary
moments on all remaining positions, (JD6) replaces the odd group's
contribution \(f/2\) by \((f-\sigma)/2\). Choose the locking and
rounding costs smaller than \(\sigma/2\). This would handle the
residual mixed-grid families that admit these moduli without a
shortness or window condition on their odd kernels. In the branch
where \(T\) has no such modulus and both square letters are wholly
nonatomic, a different joint estimate is needed unless their odd
mass is smaller than the Markov threshold.


### Conditional first moments and short retained ratios

Write \(h_3(d)=|\operatorname{Cl}(\mathbf Q(\sqrt d))[3]|\) for
squarefree \(d\), assigning value one in the trivial case. For
squarefree \(m\) and fixed \(\varepsilon,\varepsilon'>0\),
\[
 \sum_{\substack{0<|d|\le Y,\ d\ {\rm squarefree}\\m\mid d}}h_3(d)
              \ll_{\varepsilon,\varepsilon'} Y m^{-1+\varepsilon}
       \qquad(m<Y^{1/4-\varepsilon'}).
\]
There are \((h_3(d)-1)/2\) cubic fields with the same
fundamental discriminant. That discriminant is divisible by
\(m\) and has absolute value at most \(4Y\).
The cubic fields are counted by the cubic-ring divisor bound.
The trivial class contributes \(O(Y/m)\). Replacing \(d\)
by the squarefree part of \(-3d\) only changes the divisor at
\(3\) and has bounded multiplicity. The same estimate and
Cauchy–Schwarz therefore give
\[
 \#\{0<|d|\le Y:\ d\ {\rm squarefree},\ m\mid d,
                h_3(d)h_3(-3d)>n^\theta\}
        \ll_{\varepsilon,\varepsilon'}Y m^{-1+\varepsilon}n^{-\theta/2},
\]
where \(h_3(-3d)\) uses its squarefree part. This controls
sections defined by prescribing divisors of the kernel.
[BTT, Proposition 4.7](https://arxiv.org/html/2107.12819#S4.SS2)
· [Davenport–Heilbronn, Theorem 3 and its proof](https://www.math.utoronto.ca/~ila/DH2.pdf)

Conditional bounds of this kind allow concentration in a residue
fibre in an abstract dictionary. Fix \(f>\theta>0\), an integer
\(s\) sufficiently large, and \(r=\lceil s\theta/(2f)\rceil<s/2\).
For primes \(q\to\infty\), put \(n=q^{s/f}\) and take \(s\)
coordinates with labels in \(\{0,\ldots,q-1\}\). Choose independent
uniform \(b,a_1,\ldots,a_s\bmod q^r\) and set
\[
      \psi(\mathbf t)=b+\sum_{i=1}^s a_i t_i\bmod q^r,\qquad
      \mathcal D=\{\mathbf t:\psi(\mathbf t)=0\}.
\]
The images of any two distinct tuples are independent and
uniform. Some coordinate difference is a unit modulo \(q^r\),
so after fixing the other coefficients the two images determine
that coordinate's coefficient and \(b\). In a section fixing
\(k\) coordinates, the number of tuples in \(\mathcal D\) has
mean \(q^{s-k-r}\) and variance at most this mean. There are
\(\binom{s}{k}q^k\) such sections. Chebyshev's inequality with
relative tolerance \(q^{-1/4}\), summed over \(k\le\lfloor s/4\rfloor\),
therefore gives a choice for which every section has
\((1+o(1))q^{s-k-r}\) entries in \(\mathcal D\). All these limits
hold at fixed \(s\). Giving each tuple of \(\mathcal D\) the least power of three
strictly larger than \(q^r\), and giving other tuples weight one,
gives bounded first moments in every such section. Using the
same weight for both factors makes their product at least
\(n^\theta\) on \(\mathcal D\), which has size \(n^{f-rf/s+o(1)}\).

The signature factors coordinatewise in a finite unit group.
With \(X=n^{1/L}\), use \(M=q^E\) where
\(E=1+\max\{r,\lceil \nu s/(fL)\rceil\}\), and choose \(g\) of
order \(q^r\) in its principal units. The signature is
\(g^b\prod_i g^{a_it_i}\). It can also be written using any
prescribed signed odd valuation orders, by dividing the assigned
exponents by those orders modulo \(q^r\) for large \(q\).
The modulus has numerical weight at least \(\nu\) and raw
cost \(f/s\). Its identity fibre contains all of \(\mathcal D\),
whose collision energy is \(n^{2(f-rf/s)+o(1)}>n^f\) for large
\(q\). Along \(s\to\infty\), evenly spaced labels for the scaled
degrees put \(o(s)\) positions in any bounded number of shrinking
windows. The model uses assigned residue signatures and weights.
It isolates the failure of the abstract implication and leaves
the actual root congruence and class-number problem open.

For actual roots, fix a full modulus \(M\ge X^\nu\) at raw cost
\(\mu\). Omit prescribed odd positions of combined raw cost \(u\).
Fix uniform height caps \(X^{w_U+o_X(1)}\) and
\(X^{w_V+o_X(1)}\) for the retained odd products on \(U,V\).
Require \(w=w_U+w_V<\nu\) by a fixed margin. A bound on their
product separately on each row also suffices after dyadic refinement
into \(O((1+\log X)^2)\) pairs of height intervals.
Then (JD6) holds with \(\sigma=\theta/2-u\) whenever this is
positive. Fix the first odd vector at cost
\(n^{f-\theta/2+o(1)}\) by (JD5), followed by the omitted
positions of the second vector at cost \(n^{u+o(1)}\).
A collision specifies the second retained ratio modulo \(M\).
For two candidates \(P_1/Q_1,P_2/Q_2\) for this ratio, the cross
determinant \(P_1Q_2-P_2Q_1\) has absolute value \(O(X^{w+o(1)})\),
since each product uses one numerator and one denominator.
The modulus forces it to vanish. Coprimality and the exact,
coloured degree positions recover the second retained vector
from the reduced ratio. This proves the energy bound. A lock
with \(\mu<(\theta/2-u)/2\) by a fixed margin gives the saving
from (JD6). The retained odd numerical weight may be positive
and its degree measure may be diffuse.
