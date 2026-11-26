# Marginal pseudorandomness and separated tests at occurrence

At the auxiliary prime scale used to distinguish true occurrence,
an induced modular matching can have sparse marginal dictionaries
with small character sums, low additive and multiplicative moments,
and prescribed local power conditions. These properties do not
control the correlation between the three labels of a matching.
There is also a sharp obstruction for three-leg nonnegative tests
on the full product of the marginal dictionaries. Its strength
depends on the complete-letter degrees of the selected truth set.

[Support occurrence and the auxiliary-prime scale](ABC_MINIMAL_REPAIR_SUPPORT_OCCURRENCE.md#an-average-finite-field-test-equivalent-to-true-occurrence)
· [Induced additive matchings under three original radical caps](ABC_OCCURRENCE_INDUCED_ADDITIVE_MATCHINGS.md)
· [Separation rank and target multiplicity](ABC_RANDOMIZED_ADDITIVE_TEST_SEPARATION_RANK.md)

## Induced modular matchings with pseudorandom marginals

**PI1.** Fix \(5/2<\beta<3\), positive integer power orders
\(e_A,e_B,e_C\), and a finite list of integers \(h\ge2\).
Put \(k=\operatorname{lcm}(e_A,e_B,e_C)\).
For all sufficiently large odd primes \(q\), let
\(M=\lfloor q^{1/\beta}\rfloor\) and
\(H=(\mathbb F_q^\times)^k\). There are
\(m=M-O_k(M^3/q)\) rows \((a_i,b_i,c_i)\in H^3\)
whose entries within each role are distinct and for which
\[
             a_i+b_j=c_\ell
                   \quad\Longleftrightarrow\quad i=j=\ell
             \qquad\text{in }\mathbb F_q.
 \tag{PI1a}
\]
In particular, every entry meets its role's prescribed power
condition in \(\mathbb F_q^\times\).

For every one of the three marginal sets \(S\) and every
listed \(h\), the ordered additive and multiplicative
\(h\)-fold energies satisfy
\[
     E_h^+(S),\ E_h^\times(S)
               \ll_{k,h,\text{list}} M^h+\frac{M^{2h}}{q}.
 \tag{PI1b}
\]
Here \(E_h^+\) counts pairs of ordered \(h\)-tuples with equal
sum and \(E_h^\times\) counts pairs with equal product.
Both upper bounds have the order of the elementary lower
bounds from the diagonal and Cauchy–Schwarz.

Write \(e_q(x)=\exp(2\pi i x/q)\) for \(x\in\mathbb F_q\).
For every \(t\in\mathbb F_q\) and every multiplicative
character \(\chi\) of \(\mathbb F_q^\times\), one has
\[
       \left|\sum_{x\in S}\chi(x)e_q(tx)\right|
                     \ll_{k,\beta}\sqrt{M\log q},
 \tag{PI1c}
\]
except for \(t=0\) and \(\chi|_H=1\), when the sum is
\(m\) exactly. The exception includes the principal character
and every nonprincipal character trivial on \(H\).

### Construction and proof

Let \(d=[\mathbb F_q^\times:H]=\gcd(k,q-1)\).
The solutions of \(a+b=c\) with all entries in \(H\)
are invariant under simultaneous multiplication by \(H\).
There are
\[
 |H|\#\{x\in H:1-x\in H\}
   = |H|\left(\frac{q-2}{d^2}+O(\sqrt q)\right)
\]
such ordered solutions. The equality follows by setting
\(x=a/c\) and expanding the two subgroup indicators over
the \(d\) characters trivial on \(H\), extended by zero at
zero. The two principal characters contribute \(q-2\).
Every other Jacobi sum has modulus at most \(\sqrt q\),
including the degenerate character cases. This proves
nonemptiness for large \(q\) at fixed \(k\).
The action by \(H\) is transitive on the allowed values
in each role, so the three individual marginals of the
uniform law on the solutions are uniform on \(H\).
[Jacobi sum conventions and bounds](ABC_MODULAR_JACOBI_INCIDENCE.md#jacobi-sums-and-the-removal-of-bounded-order-modes)

Draw \(M\) such solutions independently, with ordered indices.
A collision of two labels in the same role occurs with probability
\(1/|H|\). A mixed equation with three distinct indices has
probability at most \(1/|H|\), by fixing its two summand rows.
When exactly two indices coincide, the mixed equation requires
a same-role collision. For example \(a_i+b_j=c_i\) forces
\(b_j=b_i\). The expected number of conflicts, counting same-role
collisions and mixed equations with three distinct indices, is
\(O_k(M^3/q)\). By Markov's inequality with a sufficiently large
constant, the conflict count is at most this order with any
prescribed fixed probability close to one.
One deletion per listed conflict suffices for (PI1a).
Because \(q\asymp M^\beta\) and \(\beta>5/2\),
the deletions cost \(o(\sqrt M)\).

For the energy bounds, first count all indexed collisions in
one marginal, before the deletions. The variables at distinct
indices on that role are independent uniform elements of \(H\).
Two \(h\)-tuples with the same index multiset contribute
\(O_h(M^h)\) collisions. Otherwise their sum difference is a
nonzero polynomial of degree one, since \(q>h\), and their
product difference is a nonzero polynomial of degree at most
\(h\). This continues to hold when some indices repeat.
A nonzero degree-\(h\) polynomial in independent uniform
\(H\)-variables vanishes with probability at most \(h/|H|\)
by the polynomial root bound on a Cartesian product.
Thus the expected indexed energy is
\(O_{k,h}(M^h+M^{2h}/q)\). A union bound with Markov's
inequality over the prescribed finite list and over the three
roles gives the bounds before deletion, with high fixed
probability. When indices and labels are made distinct, every
collision of the surviving set was already counted.
[Schwartz, polynomial root bound](https://doi.org/10.1145/322217.322225)

For (PI1c), the mean of \(\chi(x)e_q(tx)\) for uniform
\(x\in H\) is zero when \(t=0\) and \(\chi|_H\ne1\).
When \(t\ne0\), the subgroup character expansion and the
Gauss sum bound give mean \(O_k(q^{-1/2})\), uniformly in
\(\chi,t\). The total mean over \(M\) draws is
\(O_k(Mq^{-1/2})\). Apply Hoeffding's inequality separately
to the real and imaginary parts to bound the deviation by
\(C\sqrt{M\log q}\) for a sufficiently large absolute
constant \(C\), outside a probability \(O(q^{-4})\) for
each mode. A union bound over at most \(3q(q-1)\) modes
has loss \(o(1)\). The mean and the conflict deletions
both fit the bound in (PI1c). All the required events
therefore hold on one draw.
[Hoeffding, Theorem 2](https://doi.org/10.1080/01621459.1963.10500830)

### Auxiliary-prime interpretation

The prime \(q\) lies above \(M^2\), the scale where averaging
over primes in (SO1) makes false integer representatives
power-smaller than a true occurrence count of size \(M\).
In (PI1) the assertion is at a prime on residue dictionaries.
It concerns neither one fixed collection of integer lifts
on a positive proportion of primes nor a collection of lifts
satisfying the three original radical caps, the real width,
and the least-repair cutoff.

For any two of the resulting sets, independent uniform draws
have min-entropy rate tending to \(1/\beta>1/3\). They therefore
satisfy the assumptions of Lewko's Corollary 5 for the extractor
\(\mathbf 1_{\{0,\ldots,(q-1)/2\}}(xy+x^2y^2)\).
This extractor conclusion holds together with the full induced
matching in (PI1). The theorem tests independent input
samples. The arithmetic occurrence relation correlates the
samples along their common row index.
[Lewko, Corollary 5 and Lemma 6](https://arxiv.org/html/2609.27023v3#S6.SS3)

Runbo Li's multilinear prime progression result controls prime
numerators in a fixed nonzero residue class on average over
factored moduli. Its Type-II estimate assumes a rough
Siegel–Walfisz sequence. These inputs differ from test primes
dividing the varying, bounded-height errors of the actual
three-letter dictionaries. The exponent improvement in that
theorem supplies no such arithmetic lifting statement.
[Li, Theorems 1.1 and 1.4](https://arxiv.org/html/2611.29075v1#S1.SS2)

The earlier mixed-bin construction has actual primitive
integer representatives with all three radical caps and the
least-repair cutoff, at a coarser combined scale. The two
constructions address different steps of an occurrence argument.
[False mixed-bin matchings with actual capped representatives](ABC_OCCURRENCE_MIXED_BIN_FALSE_MATCHINGS.md)

## Degree-sensitive bounds for nonnegative separated tests

**PI2.** Let \(\mathcal T\ne\varnothing\) be a subset of
a finite complete product \(A\times B\times C\) and
let \(\Delta_A,\Delta_B,\Delta_C\) be its largest
one-coordinate degrees. Let \(K_\omega\in[0,1]\) on the full
product be a random acceptance kernel, with each outcome
a sum of at most \(d\ge1\) terms \(f_i(a)g_i(b)h_i(c)\)
with all factors nonnegative. If its average expected
acceptance on \(\mathcal T\) is \(p\), then
\[
 \sum_{x\notin\mathcal T}\mathbf E K_\omega(x)
  \ \ge\
       \frac{(p|\mathcal T|)^3}{d^2\Delta_A\Delta_B\Delta_C}
                             -p|\mathcal T|.
 \tag{PI2}
\]
In particular the bound uses the original complete-label
degrees when \(\mathcal T\) is an actual eligible integer
incidence relation.

For one outcome and a single term, let \(S_i\) be its sum
on \(\mathcal T\), and let \(F_i,G_i,H_i\) be the sums
of its individual factors on their respective dictionaries.
The products lie in \([0,1]\), since the kernel is bounded by
one and each term is nonnegative. Hölder gives
\[
 S_i\le \sum_{(a,b,c)\in\mathcal T}
              f_i(a)^{1/3}g_i(b)^{1/3}h_i(c)^{1/3}
     \le (\Delta_A F_i\Delta_B G_i\Delta_C H_i)^{1/3}.
\]
The full-product mass of the term is \(F_iG_iH_i\).
Sum its cubic lower bound over the terms, inserting zero
terms if needed, to obtain full mass at least
\((\sum_i S_i)^3/(d^2\Delta_A\Delta_B\Delta_C)\).
Use Jensen over the kernel outcomes and subtract the
expected true mass \(p|\mathcal T|\).

For a complete-label matching of \(J\) rows all three
degree caps equal one. If \(pJ\) is an integer divisible by
\(d\), choose \(pJ\) matched indices uniformly and partition
them into \(d\) equal groups. Accept the union of product
rectangles having the same group on all three roles.
There are exactly \(d(pJ/d)^3\) accepted cells and \(pJ\)
accepted true cells per outcome. Each true cell survives
with probability \(p\), so this attains (PI2). Unequal group
sizes differing by at most one give the analogous rounded
construction when \(pJ\ge d\). The construction here is
an extremal example for a known matching.

In the balanced exponent convention, take
\(|\mathcal T|=n^{1/3+o(1)}\), complete
dictionary exponents \(\alpha=1/3-u\), and
\(\Delta_T\le n^{u+o(1)}\) in each role, with \(p\) bounded
below. The first term in (PI2) has guaranteed exponent
at least \(3\alpha-2\log_n d\). It improves the first term
of (SR2) when \(d\) has exponent below
\(1/3-2u\) by a fixed margin after the cohort allowances.
For \(u=0\), to get a fixed normalized power below the
true-mass scale from a full-product false-mass estimate,
one therefore needs at least \(n^{1/3-o(1)}\) separated
terms. A real or arithmetic restriction on the cells counted
as noise can remove this obstruction because (PI2) prices
the mass on the full product where the decomposition is
asserted. Signed cancellation is outside its nonnegative
hypothesis.
[The complete-degree and residue cell relations (RB3–RB5)](ABC_BALANCED_DICTIONARY_RESIDUE_CELLS.md)
· [The pair-uniqueness bound (SR2) and signed flattening (SR3)](ABC_RANDOMIZED_ADDITIVE_TEST_SEPARATION_RANK.md)
