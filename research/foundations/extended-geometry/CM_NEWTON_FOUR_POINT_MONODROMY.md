# Four-point monodromy and exceptional parameters for direct-strip moments

In the fixed direct-strip examples of degrees eleven, fifteen and nineteen,
the parameter congruences give explicit certificates that the relevant
primitive Jacobian factor is not potentially CM. For each of the three
families, the same moment conclusion holds at every rational parameter
outside a finite set. The finite set is not computed here. The proof uses
complex monodromy and uniform open image for the primitive Tate module.

## Bounded-degree specializations

Let \(n\ge2\) and let \(a_0,a_1,a_\lambda,a_\infty\) be integers in
\(\{1,\ldots,n-1\}\) with sum \(2n\) and with
\(\gcd(n,a_0,a_1,a_\lambda,a_\infty)=1\). Put
\(U=\mathbb P^1_{\mathbb Q}\setminus\{0,1,\infty\}\) and let
\(C_\lambda\) be the smooth normalization of
\[
        w^n=x^{a_0}(x-1)^{a_1}(x-\lambda)^{a_\lambda}.
\]
Write \(J_\lambda^{[n]}\) for the isogeny factor of the Jacobian cut
out by the idempotent for the deck characters of exact order \(n\).

**Theorem.** For every positive integer \(D\), only finitely many
\(\lambda\in U(\overline{\mathbb Q})\) with
\([\mathbb Q(\lambda):\mathbb Q]\le D\) give a potentially CM factor
\(J_\lambda^{[n]}\).

The equation gives a connected cyclic cover over each geometric
point by the gcd hypothesis. Normalization over \(U\) gives a
smooth proper family. Locally along any branch section the Kummer
equation is a power of a uniformizer times a unit. After an étale
extraction of an \(n\)-th root of the unit, the normalized local
cover has uniformizer of ramification degree \(n/\gcd(n,a_i)\).
The four sections are disjoint, and the cover is étale elsewhere.
Thus its relative Jacobian \(J\to U\) is an abelian scheme. Over
\(\mathbb Q(\zeta_n)\) the deck generator \(\delta\) acts by
\(w\mapsto\zeta_n w\). The rational deck idempotent for exact order
\(n\) is fixed when \(\delta\) is replaced by any unit power, so
it is defined over \(\mathbb Q\). Multiplying it by an integer gives
an endomorphism whose fiberwise image represents the factor up to
isogeny.

### The rank-two period

Set \(\alpha_i=a_i/n\). On the complex family the differential
\(\eta=dx/w\) belongs to deck character \(\delta\mapsto\zeta_n^{-1}\).
At a point over the branch indexed by \(i\), its order is
\((n-a_i)/d_i-1\), where \(d_i=\gcd(n,a_i)\). At infinity this uses
\(\deg_x(w^n)=2n-a_\infty\). The order is nonnegative and strictly
less than \(n/d_i\), the order of a local parameter pulled back from
the base. Every differential in this character is a rational function
of \(x\) times \(\eta\). A pole of that function anywhere on the
base would make the differential singular on the cover. Thus the
holomorphic character space is one-dimensional. For the opposite
character the same calculation with
\(w\,dx/[x(x-1)(x-\lambda)]\) gives orders \(a_i/d_i-1\) at the
four branches and again a one-dimensional holomorphic space.
Consequently the complex character space of \(H^1\) has
\(H^{1,0}\) and \(H^{0,1}\) dimensions one and one. Its rank is two,
and the polarization gives it a monodromy-invariant
Hermitian form of signature \((1,1)\) by pairing with the conjugate
character.

For \(0<\lambda<1\) choose a continuous branch of \(w\) on the open
segment from 0 to \(\lambda\), and let \(L\) be the completed lifted
segment. The chain
\[
                    (1-\delta^{d_0})(1-\delta^{d_\lambda})L
\]
is a cycle. The stabilizer of the corresponding endpoint contains
\(\delta^{d_i}\), and the two deck operators commute. The period of
\(\eta\) on the cycle is a nonzero constant times its integral on
\(L\), because both scalars \(1-\zeta_n^{-d_i}\) are nonzero. The
integral on \(L\) itself is nonzero, as its integrand on the segment
is positive up to a constant phase. Projection onto the character
of homology detected by \(\eta\) preserves the period.

Up to the phase, substitution \(x=\lambda v\) gives
\[
 \int_L \eta
 =\lambda^{1-\alpha_0-\alpha_\lambda}
   \int_0^1 v^{-\alpha_0}(1-v)^{-\alpha_\lambda}
                    (1-\lambda v)^{-\alpha_1}\,dv.
\]
For \(|\lambda|<1\), expand the last factor by the binomial
series and integrate termwise. Dividing the beta integral in term
\(k\) by the beta integral in term zero gives
\((1-\alpha_0)_k/(2-\alpha_0-\alpha_\lambda)_k\), with
\((b)_k=b(b+1)\cdots(b+k-1)\). This yields
\[
 \int_L\eta
 = B(1-\alpha_0,1-\alpha_\lambda)\,
   \lambda^{1-\alpha_0-\alpha_\lambda}
          {}_2F_1(A,B_0;C;\lambda),
\qquad
 (A,B_0,C)=(\alpha_1,1-\alpha_0,2-\alpha_0-\alpha_\lambda),
\]
with the same phase convention. Here \(B(\cdot,\cdot)\) is the
beta function. The termwise calculation is justified by uniform
convergence of the binomial series on the interval for
\(\lambda\) in compact subsets of the unit disk and integrability
of the two endpoint weights.

The differences of the two numerator parameters \(A,B_0\) from the
two denominator parameters \(1,C\) are nonintegral. In fact
\(A,B_0,C-A=\alpha_\infty,C-B_0=1-\alpha_\lambda\) all lie in
\((0,1)\). Beukers–Heckman, Propositions 3.2–3.3, show that the
monodromy on the two-dimensional Gauss solution space is irreducible
under this condition. The prefactor \(\lambda^{1-\alpha_0-\alpha_\lambda}\)
twists it by a scalar character and preserves irreducibility.
Analytic continuations of the nonzero period thus span a space of
dimension two. Each is a period on a transported cycle in the same
rank-two character space. The period map therefore identifies this
homological representation with the twisted solution representation,
up to the convention of dualizing when transporting cohomology.
[Beukers–Heckman, §§2–3](https://math.uchicago.edu/~drinfeld/Hypergeometric_motives/Beukers-Heckman.pdf)

This complex monodromy is infinite. A finite group preserves a
positive definite Hermitian form by averaging, and in a complex
irreducible representation any invariant Hermitian form is a real
scalar multiple of that form by Schur's lemma. The form of signature
\((1,1)\) supplied by polarization rules this out. More generally,
Moonen, Proposition 4.7(i), citing Rohde, gives the special unitary
projection of connected algebraic monodromy for mixed character pairs
in families of cyclic covers of the projective line. This yields a
stronger monodromy statement than the infinitude used here. The period
calculation gives the rank-two argument explicitly for the normalized
four-point family. [Moonen, §§2 and 4](https://pure.uva.nl/ws/files/949135/83991_Moonen_DocMath_2010.pdf)

### Passage to arithmetic images

Choose a prime \(\ell\), and let \(\rho\) be the action of
\(\pi_1(U)\) on the primitive rational Tate summand, using the
lattice obtained by intersecting this summand with the full Tate
module. The intersection is stable and spans the summand over
\(\mathbb Q_\ell\). The corresponding primitive Betti space is
defined over \(\mathbb Q\), and its topological monodromy image
is infinite by the character calculation. After embedding the
algebraic closure of \(\mathbb Q\) into \(\mathbb C\), the
representations agree along the topological fundamental group by
Riemann existence and finite torsion comparison. Indeed the
analytification of \(J[\ell^k]\) is the local system of first
homology of the complex tori modulo \(\ell^k\). Passing to the
inverse limit and then applying the same rational idempotent identifies
the actions. The rational Betti matrices stay distinct over
\(\mathbb Q_\ell\), so the compact image of the geometric étale
fundamental group on the primitive summand, which contains the
closure of this topological image, is infinite.
[SGA 1, Exposé XII, §5](https://grothendiecksga.com/read/sga1/en/XII.html)

Cadoret–Tamagawa I, Theorem 5.1 in the linked version, applies to
the Tate module of the abelian scheme \(J\) over the smooth,
separated, geometrically connected curve \(U\). It gives geometric
Lie perfectness. The same holds for \(\rho\), since a subquotient
inherits this property. In particular every open subgroup of
\(G_{\mathrm{geom}}=\rho(\pi_1(U_{\overline{\mathbb Q}}))\) has
finite abelianization. If the arithmetic image \(G=\rho(\pi_1(U))\)
had an open abelian subgroup, its intersection with
\(G_{\mathrm{geom}}\) would be open and abelian, hence finite. This
would make \(G_{\mathrm{geom}}\) finite, contrary to the comparison
above. [Cadoret–Tamagawa I, Introduction and Theorem 5.1](https://webusers.imj-prg.fr/~anna.cadoret/UOI1.pdf)

At a potentially CM fiber, over a finite extension of its residue
field all elements of a commutative semisimple CM algebra of dimension
twice that of the factor are defined. The rational Tate module has
rank one over this algebra tensored with \(\mathbb Q_\ell\). This
follows from faithfulness and equality of dimensions, since every
simple component of the commutative algebra must act nontrivially.
The Galois image over the extension centralizes the algebra and
therefore lies in its commutative unit group. The specialization
image over the original residue field has an open abelian subgroup,
so it cannot be open in \(G\).

Cadoret–Tamagawa II, Theorem 1.1, applies to a geometrically Lie
perfect representation into \(\operatorname{GL}_m(\mathbb Z_\ell)\)
on a smooth, separated, geometrically connected curve over a field
finitely generated over \(\mathbb Q\). Here the field is \(\mathbb Q\)
and the chosen lattice realizes the required integral representation.
For each \(D\) it gives finitely many closed points of degree at most
\(D\) with specialization image failing to be open in \(G\). Each
closed point contributes only finitely many conjugates in
\(U(\overline{\mathbb Q})\), proving the theorem.
[Cadoret–Tamagawa II, Theorem 1.1](https://webusers.imj-prg.fr/~anna.cadoret/UOI2.pdf)

The open-image argument generalizes a finiteness argument for CM
elliptic curves in Cadoret–Tamagawa II, §4.2. A separate general
finiteness theorem of Orr–Skorobogatov gives, in any fixed dimension,
finitely many geometric isomorphism classes of CM abelian varieties
defined over fields of bounded degree. The application here concerns
parameters in the particular primitive-factor families used for moments.
[Orr–Skorobogatov, Theorem A and Theorem 2.5](https://arxiv.org/abs/1704.01647)

## Removing the parameter congruences in the direct strips

**Corollary.** In each of the three direct-strip families below, the
fixed-phase moment conclusions hold for all but finitely many rational
parameters \(c\ne0,1\).

For each row below, all but finitely many \(c\in\mathbb Q\setminus\{0,1\}\)
give a non-potentially-CM primitive factor by the theorem with \(D=1\).
The tuple lists the exponents at \(0,1,c,\infty\) in the cover
\(W^n=z^{n-1}P(z)\). The nonzero constant \(c^8\) in the first row
leaves the geometric CM property unchanged and is retained in the
prescribed angular face. The primitive coherent dimensions and the
directed paths are unchanged from the two strip theorems.

| \(n\) | Branch exponents | \(P(z)\) | Prime class | Moment order |
| ---: | --- | --- | --- | --- |
| 10 | \((9,3,3,5)\) | \(c^8(z-1)^3(z-c)^3\) | \(p\equiv7\pmod{10}\) | \((9p^2-1)/10\) |
| 14 | \((13,3,3,9)\) | \((z-1)^3(z-c)^3\) | \(p\equiv5\pmod{14}\) | \((13p^3-1)/14\) |
| 18 | \((17,2,7,10)\) | \((z-1)^2(z-c)^7\) | \(p\equiv11\pmod{18}\) | \((17p^3-1)/18\) |

For any such rational parameter, prescribe the top angular polynomial
\(z^{-1}P(z)\) on \(j+2i=n\) and require the literal support to satisfy
\((n-1)j+2ni\ge n(n-3)\). The existing strip criteria give, for
each fixed phase of this shape over a characteristic-zero field,
nonzero moments at the displayed orders on a set of rational primes
of natural density \(1/\varphi(n)\). All compatible strict lower
terms are arbitrary. The least sphere degree is \(n+1\), and the
primitive factor has \(p\)-rank zero at the certifying primes. The
finiteness proof gives no list of exceptional rational parameters.
The previous congruence certificates still give explicitly specified
subsets with no exceptions. The bounded-degree finiteness statement
also applies to algebraic parameters, where a density conclusion in
a prescribed rational prime class would require checking the residue
fields and prime-class compatibility separately.
[Degree-eleven strip](CM_NEWTON_DIRECT_STRIP_MOMENTS.md#a-degree-eleven-family)
· [Even-conductor strips](CM_NEWTON_EVEN_CONDUCTOR_STRIP_FACES.md#general-even-conductor-criterion)
