# Rank-one reduction when a prime is in the first half of a four-window

Let \(f=ax^2+bxy+cy^2+k\) with \(abck\ne0\), and
let \(M_j=\Lambda_2(f^j)\) with
\(\Lambda_2(x^i y^j)=i!j!\).

**Theorem.** Fix a prime \(p\ge7\). Suppose the four moments
at \(p,p+1,p+2,p+3\), or the four at \(p-1,p,p+1,p+2\),
all vanish. At each place over \(p\), make the coefficients
integral with at least one unit. One of the axis coefficients
is a unit. By exchanging the axes and dividing by this
unit, one can arrange that \(a=1\). In this normalization,
put \(t=b/2\) and let bars denote residue. Necessarily
\[
 v_p(k)=2/p,\qquad \overline{k^p/p^2}=-2,\qquad
 \bar c=\bar t^2,\quad \bar t\in\{0,1,-1\}.                 \tag{1}
\]
The other two positions of a prime in a four-window are
excluded by the
[prime-terminal proof](BINARY_EVEN_QUADRATIC_PRIME_TERMINAL_FOUR_WINDOW.md)
and the
[prime-penultimate proof](BINARY_EVEN_QUADRATIC_PRIME_PENULTIMATE_FOUR_WINDOW.md).
The [negative rank-one lift](BINARY_EVEN_QUADRATIC_PRIME_FRONT_NEGATIVE_RANK_LIFT.md)
excludes the residue \(\bar t=-1\) in both early-prime
positions, and the
[diagonal lift](BINARY_EVEN_QUADRATIC_PRIME_FRONT_DIAGONAL_LIFT.md)
excludes \(\bar t=1\) for the same positions. The only
possible residue is therefore \(\bar t=0\), for which the
[zero-residue conditions](BINARY_EVEN_QUADRATIC_PRIME_FRONT_ZERO_RANK_GATE.md)
give two position-specific finite-field systems. If the
corresponding two polynomials are coprime over \(\mathbb F_p\),
the four-window holds at that start. The surviving solutions
of those systems still require lifting.

The assertions can be applied over an algebraic \(p\)-adic
closure. To use them in testing existence over any
characteristic-zero field at a prescribed start, first take
an algebraic solution if one exists, by adjoining inverses
of the required nonzero expressions and using the weak
Nullstellensatz.

## An interior exclusion for a nonunit axis

Write \(Q=ax^2+bxy+cy^2\), \(H_j=\Lambda_2(Q^j)\), and
\(s=v_p(k)\). Frobenius gives \(M_p\equiv k^p\) in the
residue field, so \(s>0\). A fact used below is
\[
              H_{p-1}\equiv (b^2-4ac)^{(p-1)/2}.
                                                               \tag{2}
\]
Only the central monomial of \(Q^{p-1}\) survives the
factorial functional at this prime. Factoring the quadratic
and summing its geometric series gives (2), as in the
[terminal reduction](BINARY_EVEN_QUADRATIC_PRIME_TERMINAL_FOUR_WINDOW.md#reduction-at-the-terminal-prime).
At any order \(j\), \(H_j/(j!)^2\) is an integer-coefficient
form. Its specialization at \(Q=bxy\) is \(b^j\) and at
\(Q=a(x^2+2xy)\) is \(4^j a^j\).
[Primitive homogeneous moments](BINARY_QUADRATIC_PRIMITIVE_HOMOGENEOUS_RESULTANT.md)

Suppose both axes are nonunits, so \(b\) is a unit, or
suppose \(a=1\), \(c\) is a nonunit and \(\bar b=2\).
In both cases \(H_{p-1}\) is a unit by (2), \(H_1\) is
a unit, and \(v_p(H_{p+1})=2\) by the primitive
specializations. If \(s<1\), the only possible smallest
terms of \(M_p\) are \(k^p\) and \(p k H_{p-1}\).
The latter has valuation \(1+s<2\), the term \(H_p\)
has valuation at least two, and each other internal index
has an additional positive power of \(k\). Thus
\(M_p=0\) requires \(s=1/(p-1)\).
But at that value, in \(M_{p+1}\) the unique lowest term
comes from \(H_1 k^p\), including the unit binomial
coefficient. Its valuation \(1+s\) is less than two.
The indices \(p,p+1\) already have factorial valuation
at least two, and the indices \(2,\ldots,p-1\) each have
a factor \(p\) in the binomial coefficient and a power of
\(k\) of degree at least two. If \(s\ge1\), instead the
unique lowest term of \(M_{p+1}\) is \(H_{p+1}\), of
valuation exactly two. All other terms have valuation
at least three. This excludes both stated configurations.

## Boundary ratios across the early-prime block

After exchange and scaling we now have \(a=1\). In the
case of a nonunit \(c\), \(1-b/2\) must be a unit by the
preceding interior exclusion. If a boundary weight is
exactly zero or \(a=c\), the previously verified coefficient
planes already prove the four-window.
[Exact coupled boundary equations](BINARY_EVEN_QUADRATIC_SINGLE_FACE_STRONG_FACTORIAL.md#4-the-exact-reduction-for-general-b-and-its-limitation)

The last two indices are \(p+\ell,p+\ell+1\), with
\(\ell=1\) or \(2\). Write
\(B_j(Z)=\Lambda_1((Zx^2+k)^j)\). For unit \(Z\) and
\(1\le d\le3\), the unique lowest term in \(B_{p+d}(Z)\)
has index \(d\) when \(s<2/p\), and \(p+d\) when
\(s>2/p\). At equality these two terms alone have
valuation two. To verify this, the binomial coefficient
is a \(p\)-unit precisely for \(i\in[0,d]\) or
\(i\in[p,p+d]\). In the low unit group the smallest term
has index \(d\) and valuation \(ps\). In the high unit
group the smallest is the final one with valuation two.
At indices \(d<i<p\), the binomial contributes one
and the factorial contributes a further one when
\(i>(p-1)/2\). For \(i\le(p-1)/2\) the valuation is
\(1+(p+d-i)s\), strictly above the smaller of \(ps,2\),
including at equality. The entire range \(2d<p\) and
\(2(p+d)<3p\) holds for \(p\ge7\).

When \(c\) is a nonunit, either unique smallest term at
\(Z=1\) lies strictly below all terms at \(Z=c\); its
index is positive. Since the weight \(1-t\) is a unit,
the boundary equation rules out both strict regimes.
When \(c\) is a unit, the quotient
\(B_{p+\ell+1}(Z)/B_{p+\ell}(Z)\), after stripping its
nonzero scalar, has reduction \(Z\) in each strict regime.
It is injective on unit arguments. Distinct residues
are separated by the reduction. On each residue disk,
expand in the local coordinate about a unit point.
The monomial ratio is \(Z\), and all correction
coefficients from the higher terms in the boundary table
have positive valuation. Dividing differences by the difference of the
arguments therefore leaves a unit. The exact weighted
equations equate these two ratios at \(1,c\), which is
impossible for distinct arguments. This forces \(s=2/p\).

At \(s=2/p\), equation \(M_p=0\) forces \(\bar\Delta=0\).
Otherwise \(p k H_{p-1}\) is the unique lowest term
as in the first section, since \(1+s<ps=2\).
For the remaining restrictions write \(\kappa=\overline{k^p/p^2}\).
At the tie the coefficientwise reductions are
\[
            B_{p+d}(Z)/p^2
               \equiv (2d)! Z^d(\kappa+2Z^p),
                           \qquad 1\le d\le3.             \tag{3}
\]
The two factorial prefactors follow by Wilson's theorem.
If \(\bar c=0\), the weight of \(B_{p+d}(1)\) is a unit,
so (3) gives \(\kappa=-2\).
If \(\bar t\ne0,1\), both weights are units and
\(\bar c=\bar t^2\). The consecutive two reductions
in (3) are incompatible when \(\bar c\ne1\):
they would force both \(\kappa+2\) and
\(\kappa+2\bar c^p\) to vanish. Hence \(\bar t=-1\),
and then the weight sum reduces to \(4\), forcing
\(\kappa=-2\). Finally, when \(\bar t=1\) and
\(\kappa\ne-2\), the common factor \(\kappa+2Z^p\)
is nonzero at \(Z=1\) and cancels in the local quotient.
The ratio again reduces to a nonzero scalar times
\(Z\) on this disk, with unit linear coefficient.
This excludes equality of the exact two ratios at the
distinct lifts. The conditions (1) follow.
