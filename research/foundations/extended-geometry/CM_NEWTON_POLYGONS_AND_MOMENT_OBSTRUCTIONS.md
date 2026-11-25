# Newton-polygon bounds and sphere moment obstructions

Let \(h\) be a phase in the coordinate ring
\(K[x,y,t]/(t^2+2xy-1)\) over a characteristic-zero field. Write it in the
literal Laurent chart \(y=(1-t^2)/(2x)\). Its spherical moment is
\[
 \mathcal A(h^m)=\frac12\int_{-1}^{1}\operatorname{CT}_x(h^m)\,dt.
\]
The pure phases are those with \(\mathcal A(h^m)=0\) for every \(m>0\).
We give two exclusions in arbitrary relative rank and a transfer
criterion that uses simultaneous nonordinarity in rank two.

## A face determines a coherent Frobenius line

Suppose the highest weight \(j+\rho i\) of the literal support is
\(H=n/q>0\), with \(\gcd(n,q)=1\) and \(\rho>1\).
Suppose its initial form is specified by the Laurent polynomial
\[
 z^{-b}P(z^q), \qquad
 0<a,b<n, \quad P(0)\ne0,\quad \deg P=(a+b)/q,
\]
where the coefficient of \(z^i\) stands for the term of angle \(i\) and
height \(j=H-\rho i\). All lower-weight terms are allowed. Write the distinct
roots of \(P\) as \(\alpha_\nu\) with multiplicities \(m_\nu\), and assume
\(q m_\nu<n\) and
\(\gcd(n,n-a,n-b,qm_1,\ldots,qm_s)=1\).

For a good odd prime \(p=nk+1\), set \(m=qk\). Every term of
\(\operatorname{CT}_x(h^m)\) has height at most \(p-1\), and a term
of height \(p-1\) uses only terms on the supporting face. Even height
\(j\) contributes the denominator \(j+1\) to the integral, and odd height
contributes zero. This isolates the coefficient
\[
 p\,\mathcal A(h^{qk})\equiv
 [z^{bk}]P(z)^{qk}
 =
 [z^{p-1}](z^{n-b}P(z)^q)^k \pmod p .
\]
The congruence holds over the residue field of a degree-one prime
above \(p\) where the phase has integral coefficients and the branch
points remain distinct. The use of number-field primes also suffices
for exclusion over an arbitrary characteristic-zero field. The
coefficients and the finitely many face identities lie in a finitely
generated \(\mathbf Q\)-algebra. Localize at the leading coefficients,
root differences and omitted zero denominators, then specialize at a
maximal ideal. Its residue field is a number field. All pure moment
identities, being polynomial identities with rational coefficients, persist.

Put \(F(z)=z^{n-b}P(z)^q\) and normalize the cyclic cover \(w^n=F(z)\).
The crossing conditions give two nonzero endpoint exponents \(n-b,n-a\)
and interior exponents \(qm_\nu\), with total \(2n\). The connectedness
condition is the displayed greatest common divisor. The deck character
\(1\), defined using the action \(w\mapsto\zeta_n w\), has \(H^1(O)\)
dimension one. Indeed the finite chart has the isotypic generator \(w\)
and the infinity chart has \(w/z^2\). Both generators are integral local
generators because the reduced exponents are all less than \(n\).
The intervening Čech class is \(w/z\). Frobenius on that class sends
\(w/z\) to \(wF(z)^k/z^p\), and its scalar is the coefficient on the
right of the congruence. More generally the \(H^1(O)\) dimensions are
\[
 d_i=-1+\sum_\nu \left\{\frac{i a_\nu}{n}\right\}
\]
where this sum includes both endpoints among the branch exponents
\(a_\nu\). This is the isotypic Hasse–Witt calculation for cyclic
covers. [Bouw, §§4–5](https://doi.org/10.1023/A:1017513122376)

The primitive factor of the Jacobian is the image up to
isogeny of the rational deck idempotent for characters of order \(n\).
It carries the action of \(\mathbf Q(\zeta_n)\), with Rosati inducing
complex conjugation since the deck group preserves the Jacobian
polarization. On the primitive factor with \(s+2\) branch points its
relative rank is \(s\), because \(d_i+d_{n-i}=s\) for each unit \(i\).
This is the same normalization used in the
[primitive Image proof](../../../unpublished/primitive-balanced-image-spectrum/quintics.tex).

Thus purity makes the character-one line noninvertible under Frobenius
at every sufficiently good, split, degree-one prime. An ordinary
conjugate pair would make that line invertible. The single scalar
therefore detects nonordinarity in at least one pair. Equipping a
different number-field prime with another reduction of \(\zeta_n\)
relabels the same geometric \(w\)-line when comparing reductions of the
same model. That operation alone cannot detect a second pair in one
Newton polygon.

## One conjugate pair with multiplicities \((1,h-1)\)

Write \(h\ge2\) for the relative rank of a CM action and write \(2e\) for
the degree of the acting field. Suppose the only mixed pair has
multiplicities \((1,h-1)\), and every other pair is definite.
The part \(j=1\) of the relative norm inequality gives ordinarity after
a finite extension at all good absolute degree-one primes outside a
finite set. In fact the definite pairs contribute \(e-1\) to the sum of
ceilings of the smallest slopes. The mixed pair must therefore
contribute zero. Its height-\(h\) component with slope sum one has slopes
\(0^{h-u},(1/u)^u\) for some \(1\le u\le h\). This follows because a
positive slope \(a/d\) of multiplicity divisible by \(d\) contributes a
positive integer to the total. For \(u>1\) every slope of its polarized
dual is positive, giving a ceiling of one or more in the mixed pair.
Consequently \(u=1\), and both sides are ordinary.

This signature is the simple signature in
[Cantoral Farfán–Li–Mantovan–Pries–Tang, Assumptions 2.1 and 2.5 and Theorem 4.18](https://arxiv.org/html/2508.11174v1).
Here the arithmetic implication uses a Rosati-compatible field action
on the Jacobian factor. It therefore covers all members of the indicated
branch families without computing or requiring equality with the full
geometric endomorphism algebra.

Two examples have only one mixed pair in the primitive factor.

| Highest weight | Angular polynomial | Branch exponents | Primitive coherent multiplicities in increasing character order | Least degree of the lifted initial form |
| --- | --- | --- | --- | ---: |
| \(j+2i=4\) | \(z^{-2}P(z)\), \(P\) squarefree of degree four | \(2,1,1,1,1,2\) at order \(4\) | \(1,3\) | 6 |
| \(j+5i/3=10/3\) | \(z^{-7}P(z^3)\), \(P\) squarefree of degree three | \(3,3,3,3,8\) at order \(10\) | \(1,3,0,2\) | 8 |

In the first row the primitive factor is a fourfold with relative rank
four over \(\mathbf Q(i)\). Its displayed \((1,3)\) signature differs
from the \((2,2)\) fourfold hypothesis in
[Fité, Theorem 4](https://arxiv.org/html/2005.10185v3).
The second row gives a sixfold of relative rank three over
\(\mathbf Q(\zeta_{10})\). The lifts of the equality terms are,
respectively, scalar multiples of
\(y^2t^4,yt^4,t^4,xt^2,x^2\) and
\(y^7t,y^4t^2,yt^3,x^2\). For a negative angle \(i\), replacing the
literal leading monomial by \((-2)^{-i}y^{-i}t^{j+2i}\) adds only
strictly lower-weight Laurent terms. Nonzero initial endpoints force
the indicated least degree when the remaining normal terms have degree
at most that degree. The actual exclusion allows all strictly
lower-weight terms regardless of their degree.

The present face extraction always selects a one-dimensional coherent
eigenspace. Its pair has multiplicities \((1,h-1)\). In particular the
rank-four arithmetic signature with one \((2,2)\) pair and every other
pair definite cannot occur in the primitive factor selected by this
particular coefficient extraction.

## Transport of a vanishing line by a deck normalizer

Suppose that a Möbius map \(M\) permutes the branch set and that an
integer \(u\) prime to \(n\) satisfies
\[
 a_{M(Q)}\equiv u a_Q\pmod n
\]
at each branch point \(Q\). Then \(F(M(z))/F(z)^u\) has divisor divisible
by \(n\) on \(\mathbf P^1\). It is an \(n\)th power of a rational
function after adjoining an \(n\)th root of its constant factor.
Hence \(M\) lifts as \((z,w)\mapsto(M(z),R(z)w^u)\). This is an
automorphism since it is invertible on the base and the power \(u\)
permutes each geometric deck fiber. Its pullback sends character \(i\)
to character \(ui\) and commutes with Frobenius at the good degree-one
primes after defining the lift. Vanishing in character one transports
to every character in \(\langle u\rangle\). The count of distinct
conjugate pairs obtained in this way is
\(\big|\operatorname{im}(\langle u\rangle\to
(\mathbf Z/n)^\times/\{\pm1\})\big|\).
At relative rank two, if this count exceeds
\(\lfloor(b_{\mathrm{bal}}-1)/2\rfloor\), for
\(b_{\mathrm{bal}}>0\) balanced pairs, the partial-ordinarity bound
excludes purity. It is enough that the character-one line has
dimension one since all the transported lines then have this dimension.

For an explicit instance at the first residual signature
\(b_{\mathrm{bal}}=3\), take \(n=9\) and
\[
 P(z)=(z-1)^2(z-\omega)^3,\qquad
 \omega^2-\omega+1=0,\qquad
 F(z)=z^5 P(z).
\]
The highest weight \(j+2i=9\) with angular polynomial \(z^{-4}P(z)\)
has extreme angles \(-4,1\) and least lift degree 13. The exponent
assignments at \(0,1,\infty,\omega\) are \(5,2,8,3\).
The map \(M(z)=1/(1-z)\) cycles \(0,1,\infty\), fixes \(\omega\),
and multiplies the exponents by 4 modulo 9. The primitive dimensions
are all one for characters \(1,2,4,5,7,8\). The transported characters
\(1,4,7\) occupy all three mixed pairs. Purity would force all three
pairs to be nonordinary, and the rank-two bound allows at most one.

This example illustrates the transfer on a symmetric face. It is also
covered by the abelian-surface ordinary result. The other fixed point
of \(M\) is unramified. Composing the lift with a deck transformation
makes it fix a point over that base point, so its cube, a deck
transformation with a fixed point on an unramified fiber, is the
identity. On the primitive factor the resulting crossed product
\(\mathbf Q(\zeta_9)\rtimes\langle4\rangle\) is
\(M_3(\mathbf Q(\zeta_3))\), by the regular action on
\(\mathbf Q(\zeta_9)\) over its fixed subfield. The primitive factor
is consequently isogenous to a cube of an abelian surface.
[Sawin, Theorem 2.1 and Corollary 2.2](https://www.numdam.org/item/10.1016/j.crma.2016.01.025.pdf)
give ordinarity for that surface after finite extension.

## A single bad pair on an actual split reduction

The surviving pattern of three balanced pairs with one nonordinary
pair occurs on a reduction of the same order-nine branch family.
Replace \(\omega\) by \(\lambda=2\) and reduce at \(p=19\).
This is a good split prime, since the branch points
\(0,1,2,\infty\) are distinct and \(19\equiv1\pmod9\).
For a primitive character \(i\), let \(a_{Q,i}\) be the least residue
of \(i a_Q\) modulo 9. Normalize the finite-chart generator by
dividing \(w^i\) by \(\prod_Q (z-Q)^{\lfloor i a_Q/9\rfloor}\),
where \(Q\) runs over the finite branch points. Its ninth power has
reduced exponents \(a_{Q,i}\). As every primitive coherent dimension
is one, the same Čech calculation with \(k=2\) gives the scalar
\[
 [z^{18-2a_{0,i}}](z-1)^{2a_{1,i}}(z-2)^{2a_{2,i}} .
\]
In characters \(1,2,4,5,7,8\) the scalars modulo 19 are
\((0,7,16,7,5,0)\). The stable rank of Frobenius on the primitive coherent
cohomology is four, so the primitive sixfold has \(p\)-rank four and
the pair \((1,8)\) alone is nonordinary. This is also the pair observed
by the moment coefficient. Thus that coefficient can vanish at one
good split prime with exactly one bad pair. The example concerns a
single reduction. It makes no assertion about the density of the
pattern or about absolute degree-one primes after a chosen extension.
[Bouw, Definition 2.1 and Lemma 2.2](https://doi.org/10.1023/A:1017513122376)

## Exact finite replay

The standard-library [face checker](verify_cm_moment_cover_faces.py)
can be run with

    python3 -I verify_cm_moment_cover_faces.py

It verifies the displayed exponent assignments, coherent dimensions,
lattice lifts, order-nine orbit, and the six scalars at \(\lambda=2\).
It checks the prime coefficient
congruence by exact rational multiplication on three specified phases
with strictly lower terms added. The internal run gives isolated residues
\(6,2,18\) at \(p=13,11,19\), respectively.

It also scans connected reduced branch exponent multisets of sum \(2n\)
for \(3\le n\le30\), primitive factor dimension at most 30 and relative
rank at least three. For each multiset with one mixed pair it tries each
ordered choice of endpoints and each interior common step \(q\) coprime
to \(n\), then looks for a line with \(\rho=u/q>1\). A crossing pair
with endpoint angles \(-b,a\) supports the grid precisely when
\[
  n+ub\equiv0\pmod q,\qquad n+ub\ge2qb,\qquad n-ua\ge0.
\]
The first condition makes the heights integral. The other two impose
the endpoint lower bounds, which imply the interior lower bounds by
linearity. The endpoint degrees decrease with increasing angle.
The smallest feasible \(u\) therefore gives the least degree of a
lift for the chosen endpoints and step. This scan finds 69 branch
multisets across 20 order and rank pairs and prints one least-degree
example for each pair. Its bound concerns this face model and these
conductors. The written arguments give the exclusions for all
coefficients of the stated families.
