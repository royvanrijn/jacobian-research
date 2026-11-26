# Root projection entropy and residue occurrence

An actual dictionary of tall small-radical integers can have
independent root marginals on every small group of positions
and near-maximal growth under fixed products. Nevertheless a
small full-divisor record can be entirely determined by its
complement. An exact character-moment estimate gives some
multiscale residue control for an original box. The loss from
passing to the smaller complete dictionaries can consume that
control.

## An arithmetic dictionary from finite-field evaluations

Fix \(0<a<1/3\) as the desired complete dictionary exponent.
Use a sequence of budgets \(L=1/M\), where \(M\) tends to
infinity only after the height limit. At a fixed \(M\) put
\(n=X^L\), \(m=2M-1\), and
\[
 \beta=\frac1{3m},\qquad
 e_i=2M+i\quad(1\le i\le m),\qquad
 d=\lfloor 3am\rfloor.
\]
Take \(M\) large enough that \(1\le d<m\). The position
index \(i\) is known from its exact prime valuation \(e_i\).

Fix a prime power \(q\ge m\) once \(M\) has been fixed.
Choose distinct evaluation points \(\xi_1,\ldots,\xi_m\)
in \(\mathbb F_q\); these are fixed as \(X\) tends to infinity.
For a sufficiently small fixed \(\vartheta_M>0\), take \(m\)
disjoint pools of primes in
\[
          [(1-\vartheta_M)n^\beta,\,
                              (1-\vartheta_M/2)n^\beta].
\]
Each pool can have size \(Q_X=q^b\asymp_M n^\beta/\log X\),
where \(b=b(X)\) is an integer. For example first partition
the primes into \(m\) pools of comparable size and then take
the largest power of \(q\) supported by them. The prime
number theorem applies in this fixed relative interval at
each budget. The growing alphabet is the vector space
\(\mathbb F_q^b\); the base field and all evaluation
coordinates remain fixed in the height limit. Index each
prime pool bijectively by this alphabet.

Choose \(b\) polynomials \(f_1,\ldots,f_b\in\mathbb F_q[t]\)
of degree less than \(d\) and place the vector
\[
               (f_1(\xi_i),\ldots,f_b(\xi_i))
\]
at position \(i\). Write \(p_i\) for its indexed prime
and put \(W=\prod_i p_i^{e_i}\). The complete
integer \(W\) determines all the labels, by unique
factorization and the disjoint pools.

**RO1.** The resulting complete dictionary \(\mathcal D\)
has exactly \(Q_X^d\) values, all in \(X/2<W<X\) and all
with \(\operatorname{rad}(W)< n^{1/3}\). Its count
exponent in base \(n\) tends to \(a\) as \(L\to0\).
For any subset \(I\) of positions the full-factor
dictionary \(D_I=\prod_{i\in I}p_i^{e_i}\) has
exactly \(Q_X^{\min(|I|,d)}\) entries. With the uniform
law on \(\mathcal D\) all these projected entries are
equally likely. In particular the weak degree profile
is wholly diffuse in the outer limit and every set of
positions of total normalized radical cost at most
\(a-o(1)\) has fully saturated root marginals.

Indeed, \(\sum_i e_i=3Mm\), so \(L\beta\sum_i e_i=1\).
Take \(\vartheta_M<1/(10\sum_i e_i)\). The upper prime
bound gives \(W<X\) and the lower gives
\(W\ge (1-\vartheta_M)^{\sum e_i}X>X/2\).
The radical bound follows from \(m\beta=1/3\).
Its normalized mass measure differs by inner vanishing
errors from \(\beta\sum_i\delta_{Le_i}\), tending to
\(\tfrac13\mathrm{Unif}[2,4]\). A choice of values at
up to \(d\) distinct evaluation points gives independent
linear interpolation conditions on each \(f_j\); at \(d\)
points it determines that polynomial. This proves the
dictionary count, its uniform projection sizes, and the
limiting assertion.
[Reed and Solomon, polynomial evaluation codes](https://doi.org/10.1137/0108018)

For two subsets \(I,J\) the same
calculation gives the exact entropy composition law
\[
  H(D_I\mid D_J)
     =\bigl(\min(|I\cup J|,d)-\min(|J|,d)\bigr)\log Q_X .
                                                               \tag{RO2}
\]
In particular, if \(I\) has cost
\(t=\beta|I|\) and its complement has at least \(d\)
positions, the complete quotient \(W/D_I\) determines
\(D_I\). This is true even for packets chosen separately
on each value: the quotient identifies the indices of
the surviving positions by the disjoint prime pools and
the exact degree positions. As long as at least \(d\)
survive, interpolation recovers every omitted root.
For all such choices the conditional entropy of the
removed factor is exactly zero.

The same example has near-maximal product growth. For every
fixed \(j\ge1\), an exact \(j\)-fold product of elements of
\(\mathcal D\) determines the unordered list of the \(j\)
prime labels in every position. At each prime the product
valuation is its occurrence count times the known \(e_i\).
There are at most \((j!)^m\) ordered rearrangements of all
lists even before enforcing the same evaluation code. Thus
each product fibre has that bound, and the product image has
size at least \(|\mathcal D|^j/(j!)^m\). In particular
diagonal-scale multiplicative energy, fully saturated
small-block root projections and tall original powers
are compatible with the exact dependence in (RO2).
Distinct values agree at at most \(d-1\) positions,
since any two distinct polynomial lists differ in at
least one component. Independent values share a prime
with probability at most \(m/Q_X\), by uniformity on
a single position.

At the complete query-cell scale take \(a=u\).
The charge of every selected packet with normalized mass
\(t\) is then exactly \(\min(t,u)+o(1)\) in normalized
logarithms. Its numerical weight is at most
\(4t+o(1)\). Selecting more mass than \(u\) can give a
cheaper charge relative to the raw mass, but the record
then determines the complete query. By contrast at the
complete global dictionary scale one can take \(a=\alpha\).
These choices test different uses of the dictionary
size. The examples are complete single-letter dictionaries.
They impose neither an additive incidence condition among
three such dictionaries nor an original triple's third
radical cap.

## A residue bound at several depths

Consider a fixed-budget exact valuation box of integers
\(W=c\prod_i x_i^{e_i}\), where the coefficient \(c\) is
fixed, the \(e_i\) are fixed positive integers, and each
root position has height \(x_i\le n^{b_i}\). Zero-cost
positions are handled at height-subpower cost. The total raw
root size cap is \(\rho=\sum_i b_i\). Restrict to values
coprime to a polynomial-height modulus \(D\). Let
\(s=\log_n D\). If a selected block of positions has
total \(t>0\) with \(2t<s\) by a fixed margin, then each
unit residue class contains at most
\[
                         n^{\rho-t/2}X^\varepsilon
                                                               \tag{RO3}
\]
complete values, up to a constant allowed to depend on the
fixed budget, the positions and the margin. The selection
may depend on \(s\) after taking a finite mesh of the
normalized parameters.

Fix the roots outside the block at cost \(n^{\rho-t}\).
For a positive-cost position \(i\) inside it, take
\(h_i=\lceil t/b_i\rceil\). Then
\(h_i b_i\le2t<s\) and \(\sum_i h_i^{-1}\le1\).
Use the multiplicative characters of the unit group of
\(D\). The \(2h_i\)-moment of the sum over roots on this
position counts pairs of products whose \(e_i\)-th powers
agree modulo \(D\). The ratio of those products is one
of the \(e_i\)-th roots of unity modulo \(D\), of which
there are \(X^{o_X(1)}\) uniformly for polynomial-height
\(D\) and fixed \(e_i\). This follows from the prime-power
unit groups and the bound
\(\omega(D)=O(\log X/\log\log X)\). Each product is
strictly smaller than \(D\), so after its partner and the
root of unity are chosen it is unique. Ordered
factorizations have height-subpower multiplicity. The
moment is therefore at most \(n^{h_i b_i}X^{o_X(1)}\).
Apply Hölder to the product of character sums in the
one-class orthogonality formula. Its cost is
\(n^{t/2}X^{o_X(1)}\) and the complement freeze proves
(RO3). The argument covers in particular powers of any
fixed external prime at varying depths. All roots and
moment orders are held fixed in number at each budget.

When the individual slots are small, one can select
\(t\) proportional to \(s\) over a fixed normalized
scale range. But for an actual complete dictionary of
size \(n^\alpha\) and raw cap \(\rho\), the resulting
fractional class bound is \(n^{\rho-\alpha-t/2}\)
up to allowances. At balanced exponents
\(\alpha=\rho-u\) this is \(n^{u-t/2}\). For
\(u=\rho/2\) it supplies no decay at any scale, since
\(t\le\rho\). For a typical residue cell dictionary of
exponent \(u\le\rho/2\), the analogous fraction
\(n^{\rho-u-t/2}\) is likewise trivial. This is a
limitation of the estimate relative to those dictionaries,
even when they have arbitrary small slot granularity.
[Averaged character moments and their record cost](ABC_OVERLAPPING_MODULUS_MOMENTS.md#weighted-moments-over-a-modulus-dictionary)
· [Absolute residue bound (RB2)](ABC_BALANCED_DICTIONARY_RESIDUE_CELLS.md#regular-balanced-cohorts-and-absolute-residue-pricing)

## Scope of the projection and sieve alternatives

Shen's restricted-graph theorem works in
\(\mathbb Z/p^N\mathbb Z\) with nonconcentration across
scales for the varying \(B\) set and for the unit
coefficient set \(C\). It gives expansion for a selected
coefficient across every sufficiently dense graph on the
same two ambient sets. In a decomposition of the actual
equation such as \(A-DZ=-B\) with \(DZ\) a complete
letter, the values of \(Z\) and \(D\) vary together.
The product projections and small-block entropies of
(RO1)–(RO2) do not give a common dense input for many
coefficients. The class-mass estimate (RO3) also lacks
the relative decay required for that theorem in balanced
cells. These are separate obstacles to its use.
[Shen, Theorems 1.7–1.8](https://arxiv.org/pdf/2607.12251v1#page=4)

The November revision of Croot–Mao–Pohoata–Yip gives a
matching quadratic-intersection scale
\(\exp(c\sqrt{\log N}/\log\log N)\) at the logarithmic
prime cutoff for sets of size \(\asymp\sqrt N\) with
half-sized prime reductions. Its longer-prime-range
variant still has a subpower gain times the factor
\(|A|/\sqrt N\). The two-set variant requires both parts
near the square-root scale and a simultaneous bound
on the sum of their prime reduction sizes.
Croot–Mao–Pohoata–Sheffer–Yip's combinatorial large
sieve assumes bounded multiplicity of nonzero differences
and small reductions at every prime, and for a fixed
proportion of omitted residues its gain is again
subpower. These results do not provide the required
fixed power saving on the occupied joint residues at
the small-budget scale.
[Croot–Mao–Pohoata–Yip, Theorems 1.2–1.3, 1.5 and 1.8](https://arxiv.org/html/2607.15311v3)
· [Croot–Mao–Pohoata–Sheffer–Yip, Theorem 2.2](https://arxiv.org/html/2606.17487v1#S2.SS1)

## Injective residue registers at two fixed full steps

**RO4.** Fix numerical weights \(\nu_1,\nu_2>0\). There are
one-letter dictionaries at arbitrarily small fixed budgets \(L\)
with size \(n^{1/3-o_X(1)}\), all values in \((X/2,X)\), and radicals
below \(n^{1/3}\). Their limiting degree profiles are wholly
nonatomic. Each dictionary injects into the unit residues modulo
each of two fixed prime-power steps of numerical weights
\(\nu_1,\nu_2\). The joint register therefore has full exponent
\(1/3\) in base \(n\), and each unit progression block for either
step contains at most one dictionary value.

Take a sequence of primes \(m\ge7\) tending to infinity, and fix
\(L=m^{-2}\) before \(X\to\infty\). Put \(n=X^L\), \(J=2m-1\), and
\(b=1/(3J)\). Let \(e_j=m(2m+j)\) for \(1\le j\le J\). Choose
\(J\) pairwise disjoint prime pools in
\([(1-\tau_m)n^b,(1-\tau_m/2)n^b]\), where \(\tau_m>0\) is fixed at
this budget and \(\tau_m\sum_j e_j<1/10\). Each pool has
\(\asymp_m n^b/\log X\) primes by the prime number theorem. Form
all the values
\[
                W=\prod_{j=1}^{J}p_j^{e_j},
                \qquad p_j\text{ in pool }j .
\]
Since \(\sum_j e_j=3m^2J\), one has \(n^{b\sum_j e_j}=X\).
The prime intervals and \((1-\tau_m)^{\sum_j e_j}>9/10\)
put \(W\) in \((X/2,X)\). They also give
\(\operatorname{rad}(W)<n^{Jb}=n^{1/3}\). The choices have
distinct complete values and distinct radicals by the disjoint
pools. The normalized degree profile differs by vanishing inner
errors from \(b\sum_j\delta_{2+j/m}\), which tends to
\(\tfrac13\operatorname{Unif}[2,4]\). For any fixed integer
\(M\ge4m^2\) the actual target repair is
\[
               H_M(W)=\operatorname{rad}(W)^M
                               < n^{M/3}.
\]
Thus at the cutoff \(Y=n^{M/3}\) there are \(Y^{1/M-o_X(1)}\)
distinct repaired values in this dictionary.

Use steps \(D_1=3^{\lfloor\nu_1\log X/\log3\rfloor}\) and
\(D_2=5^{\lfloor\nu_2\log X/\log5\rfloor}\). For large \(X\)
the dictionary primes exceed five and the steps are coprime.
Write \(W=R^m\) with \(R<X^{1/m}\) and choose \(m\) large enough
that \(1/m<\min(\nu_1,\nu_2)\). Then \(R<D_1,D_2\) at
sufficiently large heights. The \(m\)-th power map is a
bijection on each unit group since \(m\) is prime to both
Euler totients. Equality of two \(W\) residues modulo either
step forces equality of their \(R\) residues and then of their
positive roots. Each residue map is therefore injective.
The construction concerns one complete role. The missing
estimate couples the occupied registers to the other two roles
of the original additive equation.
[Related exact-power dictionaries with diffuse degrees](ABC_OCCURRENCE_INDUCED_ADDITIVE_MATCHINGS.md#a-primitive-two-cap-induced-matching-with-diffuse-degrees)
