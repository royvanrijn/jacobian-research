# Even-conductor direct strip faces in degrees fifteen and nineteen

An even cyclic conductor with one balanced primitive pair gives an
odd source character for the directed moment test. The source line
has a specified generator, and the same strip works for every
compatible lower term. Besides the order-ten construction, this
gives fixed phase families at conductors 14 and 18, with moments
detected on rational prime sets of density \(1/6\). These new
conductors and degrees arise from different primitive Jacobian
factors.

## General even-conductor criterion

Let \(n\ge4\) be even. Suppose \(P(z)\in\mathbb Q[z]\) has exactly
two distinct nonzero geometric roots, with positive multiplicities
and total degree \(a+1\), where \(1\le a\le n/2\).
Consider the normalized curve
\[
                         W^n=z^{n-1}P(z).
\]
It is connected since one branch exponent is \(n-1\).
The exponents at zero and infinity are \(n-1,n-a\); the other
two are the two multiplicities of \(P\). Thus every oriented
primitive character has relative rank two, and character one
and its conjugate \(n-1\) have coherent dimension one by the
fractional-part formula.

Assume the primitive factor is not potentially CM and the pair
\(\{1,n-1\}\) is its only balanced conjugate pair. Suppose a unit
\(u\bmod n\) has some \(r\ge2\) with \(u^r\equiv-1\bmod n\),
with \(r\) the least such exponent, and that the intermediate
characters \(-u,-u^2,\ldots,-u^{r-1}\) all have nonzero coherent
dimension. Every one of those dimensions is then two.
Over \(\mathbb Q\) the primitive factor carries the geometric
action by \(E=\mathbb Q(\zeta_n)\) as a Galois-stable algebra,
since the primitive rational deck idempotent is invariant.
Rosati on \(E\) is conjugation. The
[one-balanced-pair theorem](CM_NEWTON_SEMILINEAR_PROJECTIVE_BOUNDS.md#averaged-hodge-types-and-a-single-balanced-pair)
and the [selected-path lemma](CM_NEWTON_MU_ORDINARY_SELECTED_PATHS.md#a-segment-with-nonzero-coherent-dimensions)
give a nonzero map of lines from character \(n-1\) to character
one on a density-one subset of good primes in the class
\(p\equiv u\bmod n\). The return on character one is zero:
its next character \(u\) has dimension zero by duality with
\(-u\). If \(u\) generates the primitive character group, every
primitive character lies on this same Frobenius orbit. The
primitive coherent Frobenius is then nilpotent and the primitive
factor has \(p\)-rank zero.

The source class used in the moment is itself a generator.
If the two nonzero roots are \(\alpha_\nu\) with multiplicities
\(m_\nu\), the normalized finite generator of character \(n-1\)
is
\[
 e_{n-1}=\frac{W^{n-1}}{z^{n-2}
                       \prod_{\nu=1}^2(z-\alpha_\nu)^{m_\nu-1}} .
\]
This uses \(\lfloor(n-1)e/n\rfloor=e-1\) at each positive
reduced finite branch exponent \(e\). Hence
\[
 [W^{n-1}/z^{n-1}]
  =\left(\prod_{\nu=1}^2(-\alpha_\nu)^{m_\nu-1}\right)[e_{n-1}/z].
\]
The remaining terms are regular on the finite chart. The scalar
in parentheses is nonzero away from finitely many primes.

**Theorem (strip conclusion).** Fix a sphere phase over any characteristic-zero
field with largest literal Laurent weight \(j+2i=n\), angular
polynomial \(z^{-1}P(z)\) on that line, and all literal support
satisfying
\[
                  (n-1)j+2ni\ge n(n-3).
\]
All strictly lower-weight terms satisfying the strip are arbitrary.
Its least sphere degree is \(n+1\), and
\[
                 \mathcal A\bigl(f^{((n-1)p^r-1)/n}\bigr)\ne0
\]
on a set of rational primes in \(p\equiv u\bmod n\) of natural
density \(1/\varphi(n)\). The phase is fixed throughout this prime
set. The directed coherent map certifies the moments even though
the selected return is zero.

Indeed this is the [strip criterion](CM_NEWTON_DIRECT_STRIP_MOMENTS.md#a-strip-criterion)
with \(s=n-1\), \(B=n(n-3)/(n-1)\), and highest weight \(n\).
Equivalently for \(m=((n-1)p^r-1)/n\), any angular-constant
height is at least
\[
                 Bm=(n-3)p^r-\frac{n-3}{n-1}
\]
and at most \((n-1)p^r-1\). The denominator \((n-2)p^r\)
has odd height and integrates to zero. Hence multiplying by
\((n-1)p^r\) and reducing leaves \([z^m]P(z)^m\), the selected
Frobenius entry just proved nonzero. Primes dividing \(n-1\)
and any bad denominators are excluded. For coefficients over an
arbitrary characteristic-zero field, specialize the finite algebra
of coefficients to a number field, preserving the leading face.
There are only finitely many additional primes of nonintegrality
after specialization. Nonvanishing of the specialized moment
implies nonvanishing of the original one.

For the exact degree, the negative equality angle \(-1\) lifts
by \(-2y t^n=x^{-1}t^{n+2}-x^{-1}t^n\). Both terms are in
the strip. Nonnegative equality angles \(0,\ldots,a\) have
nonnegative heights because \(a\le n/2\), and lie in the strip.
At any negative angle \(-v\), the spread of allowed \(t\)-heights
is at most \(2(n-v)/(n-1)\), but regularity requires a spread
at least \(2v\) from the factor \((1-t^2)^v\).
Thus \(v\le1\). In the normal sphere basis the negative
angle can have degree at most \(n+1\), while nonnegative
angles have degree at most \(n\). The nonzero extremal
\(x^{-1}t^{n+2}\) requires a term \(y t^n\), so the degree
is exactly \(n+1\).

## Two new conductors

**Corollary.** The following data satisfy the hypotheses. The polynomial is
\(P(z)=(z-1)^{m_1}(z-c)^{m_2}\) in each row.

| Order \(n\) | \((m_1,m_2)\) | \(c\) modulo the two indicated primes | Prime class | Coherent path | Moment order |
| ---: | --- | --- | --- | --- | --- |
| 14 | \((3,3)\) | \(3\) at 43 and 127 | \(5\bmod14\) | \(13\to9\to3\to1\) | \((13p^3-1)/14\) |
| 18 | \((2,7)\) | \(3\) at 19 and 109 | \(11\bmod18\) | \(17\to7\to5\to1\) | \((17p^3-1)/18\) |

Take rational \(c\) integral at the primes indicated. The
primitive coherent dimensions in increasing character order
are \(1,2,0,2,0,1\) at 14 and \(1,2,2,0,0,1\) at 18.
In both rows the path has dimensions \(1,2,2,1\).
Thus the two families give nonzero moments of fixed phases
of least sphere degrees fifteen and nineteen, each on a set of
natural density \(1/6\). The residue class in each row generates
the six primitive characters, so the primitive factors have
\(p\)-rank zero at the certifying primes. The strip imposes no restrictions
on lower coefficients beyond its support inequalities.

Here is the finite-field certificate that the primitive factors
are not potentially CM. Let \(\Delta=a_p^2-4\delta_p\) be
the relative Frobenius discriminant of the primitive
rank-two action over \(E=\mathbb Q(\zeta_n)\). The exact
character-sum calculation over \(\mathbb F_p\) and
\(\mathbb F_{p^2}\) gives:

| \(n\) | \(p\) | \(N_E(\Delta)\) | Auxiliary prime for the distinct-power test |
| ---: | ---: | --- | ---: |
| 14 | 43 | \(43^4\,1583^2\) | 109 |
| 14 | 127 | \(13^2\,71^2\,127^4\) | 109 |
| 18 | 19 | \(3^4\,19^4\,71^2\) | 103 |
| 18 | 109 | \(19^2\,109^4\,163^2\) | 103 |

In each row the norm of the relative polynomial is irreducible
of degree twelve over \(\mathbb Q\). No eigenvalue ratio is
a root of unity: such a ratio would have degree at most 144,
and its order would be killed by
\(\operatorname{lcm}\{j:\varphi(j)\le144\}\).
The checker verifies that the polynomial of eigenvalue powers
at exactly this exponent is squarefree at the indicated
auxiliary prime. The search for \(j\) is bounded by \(41472\)
using \(\varphi(j)\ge\sqrt{j/2}\).
Thus all the finite reductions used here are absolutely simple,
with geometric endomorphism algebra the degree-twelve Frobenius
field, by Tate's theorem.

At \(n=14,p=43\) two of the split places above 1583 divide
\(\Delta\), corresponding to the values 261 and 279 of
\(\zeta_{14}\) modulo 1583. The norm valuation is two,
so both valuations are one. At \(p=127\) the relative
discriminant is a unit at 1583. At \(n=18,p=19\), 71 has
residue degree two in \(\mathbb Q(\zeta_{18})\), so the norm
valuation two forces a valuation-one discriminant place.
At \(p=109\), 71 does not divide the norm.
In each conductor, the two geometric reduction fields therefore
differ in ramification. A characteristic-zero geometric CM
algebra of dimension twelve could not inject into both.
The same certificates apply to every parameter in the table,
since the two required reductions are preserved.
[Exact Sage replay](verify_cm_newton_even_conductor_strip_faces.py)
· [Tate, Theorem 1](https://numdam.org/item/SB_1968-1969__11__95_0/)

The parameter congruences allow infinitely many rational cross-ratios
in each conductor. The arithmetic statement applies to each fixed
parameter and each compatible fixed phase. It does not require a
uniform exceptional set across parameters.

The internal SageMath replay checks all four arithmetic records,
the two full primitive Frobenius orbits, and the selected coefficient
at \(p=5,m=116\) for order 14 and \(p=29,m=23034\) for order 18.
From this directory use a Python interpreter with SageMath available.

    python verify_cm_newton_even_conductor_strip_faces.py
