# Exact tame quotient separation under arbitrary factorial tails in repeated blocks

Fix an ordinary two-monomial exponent profile \(C\). For an adjacent
pair of orders write \(n,N=n+1\), \(n\ge2\). Let \(c_{j,k}\) be the
\(p\)-adic valuation of the ordinary coefficient at index \(k\)
after removing the constant coefficient, where \(j=n,N\).
Fix a prime \(p>N\). The binomial coefficients are then all units.

There is a finite hyperplane test on the two vectors \(c_j\).
It has the following consequences.

**Theorem.** Suppose the two vectors pass the test modulo a
positive integer \(L\). Append any number of ordinary exponent
rows, with every appended row type having multiplicity divisible
by \(L\). The resulting adjacent moments have disjoint root
arguments. The appended gaps, signs, heights, and total curvature
are arbitrary. The hypothesis gives exact separation near every
positive exceptional Gaussian strength of a marked core as well
as at all other strengths.

Let \(J\ge4\) be an integer at least the
[absolute ordinary terminal cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
One can choose \(p,L\), depending only on \(J\), so that every
ordinary core has an extension by at most \(2\sum_{j=2}^{J-1}(j-1)\)
unit-gap rows whose valuation vectors pass simultaneously for all
\(3\le N<J\). The added rows can have arbitrarily small total
curvature and maximum relative gap, preserve the signed total gap
exactly, and make the physical coefficient multipliers after
constant normalization arbitrarily close to one at any prescribed
finite list of orders. In particular any subsequent tail in
multiplicity blocks of \(L\) gives ordinary SFC(2) in every degree,
uniformly over its counts and heights.

The certificate concerns the actual integer coefficients. It uses
neither separation of a Gaussian reference nor a lower bound on
how far an ordinary tail is from an exceptional strength.

## Finite valuation test

For a real vector \(v=(v_0,\dots,v_j)\) and an edge candidate
\(I=[l,r]\), \(0\le l<r\le j\), define
\[
 s_I(v)=\frac{v_r-v_l}{r-l},\qquad
 d_{I,k}(v)=v_k-v_l-(k-l)s_I(v)\quad(k\notin\{l,r\}).
\]
The forms below use the vector from the same order as the edge.
Take \(D=\operatorname{lcm}(1,\ldots,N)\). The test requires each
of the following linear forms, after multiplication by \(D\),
to be nonzero modulo \(L\).

1. In either order, every \(d_{I,k}\), and every difference
   \(d_{I,k}-d_{I,k'}\) with \(k\ne k'\).
2. Across the two orders, every difference \(d_{I,k}-d_{I',k'}\).
3. Across the two orders, every difference
   \[
      \bigl(s_{I_2}-s_{I_1}\bigr)
             -\bigl(s_{I'_2}-s_{I'_1}\bigr),
   \]
   where each pair of edges is ordered left to right with
   disjoint interiors. The two edges of a pair may share a
   vertex.

All forms have integer coefficients after clearing denominators.
Each is invariant under an independent affine addition to each
vector. Every form is nonzero on the space of unrestricted
vectors modulo those affine additions. For a difference of
slope differences, each within-order difference is itself
nonzero, as can be seen at the first endpoint. For differences
of two depths in a single order, their coefficients at the
two nonendpoint indices already show nonvanishing. The other
cases each contain a nonzero coefficient at a nonendpoint
index from one order.

In a tail whose row multiplicities are multiples of \(L\),
each stripped valuation changes by a multiple of \(L\).
Consequently all three tests remain true for the exact
coefficients of the extended profile. In particular the
corresponding rational forms are nonzero there.

## What one tame prime detects

Use the extended profile's actual valuations \(v\).
Every edge of its lower Newton hull has exactly two sampled
points on its supporting line by test 1. Let \(I=[l,r]\) be
one such edge. Scale the polynomial variable by an element
of valuation \(-s_I(v)\), and divide by a nonzero scalar so
that both face coefficients are units. The scaled polynomial
\(F\) has integral coefficients, reductions zero at all
off-face indices, and unit face coefficients at \(l,r\).
Its reduction has \(r-l\) distinct nonzero roots, since
\(p>N\). These reduction roots lift uniquely. All roots
with valuation \(-s_I(v)\) are accounted for by these lifts.

Let \(\alpha,\beta\) be different unit roots of \(F\)
from this edge. Their reductions differ by a nontrivial
root of unity of order dividing \(r-l\). Let \(\omega\) be
the unique root of unity of that order with
\(\bar\omega=\bar\beta/\bar\alpha\). Then either
\(\beta/\alpha=\omega\), or
\[
 v_p(\beta/\alpha-\omega)
    =\min_{\substack{k\notin\{l,r\}\\\omega^{k-l}\ne 1}}
                   d_{I,k}(v).                                      \tag{1}
\]
An empty minimum is interpreted as \(+\infty\) and occurs
precisely in the first case. In fact, if \(f_k=[z^k]F\),
\[
 F(\omega\alpha)-\omega^l F(\alpha)
     =\sum_k f_k \omega^l(\omega^{k-l}-1)\alpha^k .
\]
The face terms and every invariant off-face term vanish.
Every other \(\omega^{k-l}-1\) is a unit, because its order
is prime to \(p\). The remaining terms have distinct
valuations by test 1, so their smallest valuation is the
valuation of the sum. The divided difference between
\(F(\omega\alpha)\) and \(F(\beta)=0\) is a unit, as both
arguments reduce to the same simple nonzero root. This
proves (1), including the case of an empty sum. For the
claim about residue roots one may work in a finite extension
of \(\mathbb Q_p\) containing the scaling elements, division
scalars, roots of unity, and polynomial roots, and then
enlarge its residue field to split the face binomials.

A quotient of roots from distinct edges has nonzero valuation.
Its absolute valuation is the difference of the corresponding
increasing slopes. Test 3 gives disjoint nonzero absolute
quotient valuations across the two orders. A quotient from
a single edge has valuation zero. If such quotients in the
two orders were equal, they would reduce to the same
nontrivial root of unity of order \(h\) dividing both chosen
edge lengths. When both minima in (1) are finite, their
values differ by test 2. If exactly one minimum is infinite
the quotients also differ. If both minima are infinite,
every index in each order is congruent to the left endpoint
of its chosen edge modulo \(h\). Each array has every
coefficient nonzero, in particular at the consecutive
indices \(0,1\), so this is impossible for \(h>1\).

We have proved disjointness of all off-diagonal root
quotients of the actual polynomials. To check arguments
over \(\mathbb C\), a shared nonreal argument would give
equal quotients of each root by its complex conjugate.
These algebraic equalities would hold under an embedding
into \(\overline{\mathbb Q}_p\), contradicting the
valuation test. The even moment is strictly positive on
the real axis for distinct monomials. Thus no shared
argument is real either.

## Simultaneous ordinary realizations at one prime

Take \(p>J^2\) and use the distinct rational numbers \(k/j\)
with \(2\le j<J,\ 1\le k<j\). For any such fraction \(f\),
an increasing row \((A,A+1)\) has stripped successive
factorial ratios \(jA+k\). Impose \(A\equiv-f\pmod p\).
No other distinct fraction on these grids has that residue,
because a nonzero difference of two grid fractions has
a numerator of absolute value less than \(J^2\) before
reduction. The terminal ratio \(k=j\) is also unmarked.
For every desired positive integer \(e\), one can prescribe
\(v_p(A+f)=e\) via a residue class modulo \(p^{e+1}\).
Hence one row with that residue adds valuation \(e\)
exactly at the successive-ratio coordinate \(f\) in every
grid where \(f\) occurs. If \(e=0\) omit the row.

In adjacent orders the open fractional grids are disjoint,
since \(n,N\) are coprime. Their internal successive-ratio
ordinates, with the terminal ratio fixed, are independent
coordinates for the two vectors modulo affine additions.
Specifically, setting the constant ordinate to zero and
the last successive-ratio valuation to zero fixes the
two affine degrees of freedom for each order. Thus each
of the finitely many cleared forms for all orders remains a
nonzero integral linear form in the distinct fractional-grid
coordinates. Pick \(L\) prime larger than the number of
forms and large enough that every form still has a nonzero
coefficient modulo \(L\). Adding the fixed core contributes
only a constant to each constraint. Each forbidden affine
hyperplane over \(\mathbb F_L\) occupies a fraction \(1/L\)
of the residue space. Some residue vector avoids all of
them. The residue search is small-dimensional for fixed
\(J\), and existence uses just the union bound.

Realize its coordinates by the increasing rows above,
using representatives \(0,\ldots,L-1\). For each realized
row use a decreasing unit-gap companion \((B+1,B)\)
where \(B\) is a multiple of \(p\) with \(|B-A|<p\).
All its stripped factorial factors in orders \(j<J\) are
\(p\)-units, because they are among \(jB+1,\ldots,jB+j\).
The signed total gap of each pair is zero. Its stripped
coefficient multiplier, in the physical variable without
any argument rescaling, is exactly
\[
        \prod_{k'=1}^{k}
                  \frac{jA+k'}{jB+j-k+k'}.
\]
It tends uniformly to one on each fixed finite list of
orders as \(A\) grows in its admissible residue class.
The pair's curvature and relative gaps go to zero at
the same time. There are at most as many pairs as
fractional-grid coordinates. The finite valuation tests,
positivity at terminal two, and the absolute angular
cutoff prove the all-degree assertion.

## Arbitrarily many repetitions of finitely many fixed blocks

Fix a core and any \(t\) prescribed ordinary row multisets
\(B_1,\ldots,B_t\), each with any finite number of rows.
An extension of the core by the same bound on the number of
unit-gap rows can certify all repetitions
\[
               C^{\rm ext}+\sum_{i=1}^t r_i B_i,\qquad
               r_i\in\mathbb Z_{\ge0},
\]
in every degree at once. The added rows still have zero
total signed gap, arbitrarily small curvature, and arbitrarily
small change to the physical coefficient ratios on finite
order grids. Arbitrary additional rows repeated in multiples
of \(L\) are allowed. The modulus \(L\) in this assertion
can be chosen to depend only on \(J,t\).

Indeed, use \(t+1\) distinct primes \(p_a>J^2\) and impose
independent marker residues modulo \(L\) in the fractional
grid at each prime. Let \(m\) be the total number of cleared
forms for all adjacent orders below \(J\), including the
order labels. Pick \(L>m^{t+1}\) prime with each form still
nonzero modulo \(L\). If every prime fails for a repetition
vector \(r\) modulo \(L\), choose a failing form \(F_a\)
at each prime. For this fixed choice, the \(t+1\) vanishing
conditions form a system of affine equations in the \(t\)
unknown residues of \(r\). Its matrix is fixed by the
valuation arrays of the \(B_i\). As the marker vectors
vary independently and uniformly at the different primes,
the vector of constant terms is uniform in
\(\mathbb F_L^{t+1}\). The matrix has rank at most \(t\),
so the probability the system is consistent is at most
\(1/L\). The union bound over all \(m^{t+1}\) choices
gives a set of marker vectors for which some prime passes
all forms at every \(r\). The proof can be made deterministic
by solving these finite linear systems for candidate marker
vectors.

For a given grid fraction use one increasing unit-gap row
to realize all its chosen depths at once. At a prime with
prescribed depth zero impose \(A\equiv0\); at each other
prime impose the earlier residue for its desired depth.
The Chinese remainder theorem gives arbitrarily large
simultaneous solutions \(A\). Choose its decreasing
companion with \(B\) a nearby multiple of the product
of the \(p_a\). This preserves the row bound and all the
physical and curvature claims. The full certificate is
again a statement about the original factorial coefficients
of each indicated extension.

The modular test is also applicable directly to an original
profile, with no added rows. These transfers rely on
congruence control on both full coefficient arrays.
An unrestricted tail with independently chosen exponents
and multiplicities does not come with such control.

## Why finitely many fixed primes cannot force this test

At one chosen adjacent pair, fix any finite list of primes
\(p_a>N^2\). From any ordinary profile one can append
signed-gap-zero pairs with arbitrarily small curvature and
physical coefficient multipliers arbitrarily close to one
in both orders. The exact extended valuation
vectors at all the chosen primes can be made affine in both
orders. In particular the modular tests at those primes all
fail for the same profile.

To see this, at a prime \(p_a\) let \(u_{j,k}\) be the
profile's valuations of successive normalized coefficient
ratios. Choose a nonnegative integer \(E\) large enough that
\(u_{j,j}+E\ge u_{j,k}\) in each order for every \(k<j\).
Mark fraction \(1\) at depth \(E\), increasing both last
ratios by \(E\). The remaining marks are independent on
the two open grids. At fraction \(k/j\) choose depth
\(u_{j,j}+E-u_{j,k}\). Every step valuation in each order
is now the same. Coordinatewise Chinese remaindering over
the chosen prime list and addition of the \(p_a\)-unit
decreasing companions preserve exactly the required
depths and the zero signed added gap. Their arbitrarily
large admissible heights give the small curvature and
physical coefficient changes from the previous section.

At a further prime \(p_*>J^2\) the same construction can
prescribe simultaneous passing residue patterns at all
orders below \(J\). Use the union of the fractional grids
and give its unused coordinates depth zero at a prime
where they are not wanted. In the selected adjacent pair
they cannot affect a prescribed factor because its
fractions have distinct residues modulo \(p_a\).
Thus the flat examples can themselves be taken all-degree
coprime. The construction forces failure of these valuation tests
on actual all-degree coprime profiles. It holds even for extensions of
already diffuse tails arbitrarily close to an exceptional
strength, since the added pairs can have arbitrarily
small curvature, zero signed gap, and zero signed
midpoint bias in the limit.
