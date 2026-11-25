# Fixed leading faces for nonsplit directed Frobenius moments

The directed coherent Frobenius obstruction can test a fixed
unmarked phase on a positive-density set of nonsplit primes.
The phase needs no radial term depending on the prime. A lower
linear bound on the literal Laurent support removes the smaller
interval poles, and the parity of interval integration removes
the intervening pole.

## A strip criterion

Use the literal Laurent chart \(y=(1-t^2)/(2x)\) on the sphere.
Suppose the highest \(j+\rho i\) weight of a phase is \(n\),
with angular polynomial \(z^{-b}P(z)\) there, and suppose all
literal monomials \(x^it^j\) also satisfy
\[
                   j+\omega i\ge B .
\]
Assume \(s>1\) is odd and \(B\ge n(s-2)/s\). For an odd prime
\(p\nmid s\), take \(r\ge1\) for which \(s p^r\equiv1\bmod n\),
and let \(m=(s p^r-1)/n\). An angular-constant product has
height \(j\) in the interval \([Bm,nm]\), because the angular
exponents sum to zero in both support bounds. In particular
\[
                     (s-2)p^r-1<j\le s p^r-1 .
\]
An interval denominator \(j+1\) divisible by \(p^r\) is either
\((s-1)p^r\) or \(s p^r\). The former has odd height and
integrates to zero. Hence for integral coefficients at the prime,
\[
 s p^r\mathcal A(f^m)
   \equiv [z^{bm}]P(z)^m\pmod p .
\]
Only the leading weight can participate at the surviving height.
There is no weight congruence on the smaller terms. This complements
the [support-congruence method](CM_NEWTON_SELECTED_EIGENLINES.md#a-moment-detecting-a-transition-rather-than-a-return).

For the connected Kummer cover \(W^n=z^{n-b}P(z)\) with reduced
nonzero branch exponents and total branch exponent \(2n\),
the right side is the selected entry of \(F^r\) from character
\(s\bmod n\) to character one. Indeed the overlap class of
\(W^s/z^s\) is sent to \(W(z^{n-b}P(z))^m/z^{s p^r}\);
the coefficient at the character-one line generator \(W/z\)
is precisely the one in the congruence.
[Bouw, §§4–5](https://doi.org/10.1023/A:1017513122376)
· [Lin–Mantovan–Singhal, §2.4](https://arxiv.org/html/2406.09632v1#S2.SS4)

## A degree-eleven family

**Theorem.** Fix \(c\in\mathbb Q\), integral at 11 and 31 with
\(c\equiv3\) at both primes. Let \(f\) be any phase over a
characteristic-zero field whose literal support satisfies
\[
                    j+2i\le10,\qquad 9j+20i\ge70 .
\]
Prescribe its angular polynomial on \(j+2i=10\) to be
\[
              z^{-1}P(z),\qquad
              P(z)=c^8(z-1)^3(z-c)^3 .
\]
All strict lower-weight terms in the strip are arbitrary.
For a set of rational primes \(p\equiv7\bmod10\) of natural
density \(1/4\),
\[
                     \mathcal A\bigl(f^{(9p^2-1)/10}\bigr)\ne0 .
\]
At the certifying primes the primitive order-ten factor has
\(p\)-rank zero. The conclusion is for each fixed parameter and phase.

**Proof.** First the face has regular sphere lifts. Its only
negative extreme angle is \(-1\), with leading height 12.
The normalized lift \(-2y t^{10}\) has the corresponding
literal terms \(x^{-1}t^{12}-x^{-1}t^{10}\).
Both satisfy \(9j+20i\ge70\). All nonnegative equality
angles \(0,\ldots,5\) also satisfy it. The degree assertion
is proved for the scaled strip below.

Consider the fixed curve
\[
             D_c:\quad w^{10}=u^5(u-1)^3(u-c)^3 .
\]
Its smooth normalization has genus seven, with branch
exponents \(5,3,3,9\). The primitive order-ten Jacobian factor
has dimension four, is defined up to isogeny over \(\mathbb Q\),
and carries \(E=\mathbb Q(\zeta_{10})\). Indeed the rational
deck idempotent for primitive characters is Galois invariant.
Rosati inverts the polarized deck transformation, and Galois
preserves the \(E\)-action as an algebra. By the fractional-part
rule the primitive coherent multiplicities in characters
\(1,3,7,9\) are \(1,2,0,1\), so there is a single balanced pair.

The primitive factor is not potentially CM. Its good reductions
at 11 and 31 are the same for every permitted \(c\), since
the branch parameters reduce to 3. Write \(\zeta=\zeta_{10}\).
Order-ten character sums on these residue fields and their
quadratic extensions give the following trace and determinant
of the relative degree-two Frobenius polynomial:
\[
\begin{array}{c|l|l|l}
p&a&\delta&N_E(a^2-4\delta)\\\hline
11&\zeta^3+\zeta^2-\zeta&
  -9\zeta^3+3\zeta^2-7\zeta-3&11^2\,149^2\\
31&-\zeta^3+2\zeta-3&
  -17\zeta^3+20\zeta^2-26\zeta+38&31^2\,389^2
\end{array}
\]
The norms of the relative polynomials are irreducible
over \(\mathbb Q\), giving degree-eight Frobenius fields
containing \(E\). At 11 the field ramifies above 149:
the residue degree of 149 in \(E\) is two and the displayed
norm valuation of the relative discriminant is two, so
there is a valuation-one place. At 31 the corresponding
field is unramified above 149 since its relative discriminant
is a unit there.

The reductions are absolutely simple. If a ratio of two of
their eigenvalues were a root of unity of order \(j\), the
eigenvalues together generate a field of degree at most 64,
so \(\varphi(j)\le64\) and \(j\le8192\), using
\(\varphi(j)\ge\sqrt{j/2}\). The characteristic polynomial of
the eigenvalue powers with exponent
\(N=\operatorname{lcm}\{j:\varphi(j)\le64\}\) is squarefree
modulo 137 at 11 and modulo 103 at 31. Hence the eigenvalues
remain distinct under every power, and irreducibility persists
over every finite residue extension. Tate's theorem identifies
the geometric endomorphism algebras of the reductions with
the degree-eight fields just found. A geometric CM algebra of
dimension eight in characteristic zero would inject into both
fields by good specialization, impossible by their ramification.
[Exact Sage replay](verify_cm_newton_direct_strip_moments.py)
· [Tate, Theorem 1](https://numdam.org/item/SB_1968-1969__11__95_0/)

The [one-balanced-pair theorem](CM_NEWTON_SEMILINEAR_PROJECTIVE_BOUNDS.md#averaged-hodge-types-and-a-single-balanced-pair)
now gives \(\mu\)-ordinary reduction relative to \(E\) at a
density-one set of rational primes. Abel–Jacobi identifies the
curve's coherent cohomology with that of its Jacobian, compatibly
with deck actions and Frobenius; the primitive idempotent retains
the indicated eigenspaces. In the class \(p\equiv7\bmod10\)
the coherent segment \(9\to3\to1\) has dimensions \(1,2,1\).
The full primitive orbit passes through character seven of coherent
dimension zero, so the coherent Frobenius is nilpotent on the primitive
factor and the stable-rank characterization gives \(p\)-rank zero. The
[selected-path lemma](CM_NEWTON_MU_ORDINARY_SELECTED_PATHS.md#a-segment-with-nonzero-coherent-dimensions)
makes the directed map of lines nonzero.

Use the rational change \(u=c/z\), \(W=z^2w\). The equation
becomes \(W^{10}=z^9P(z)\), precisely the cover of the leading
face. Its normalized finite generator at character nine is
\[
                    e_9=\frac{W^9}{z^8(z-1)^2(z-c)^2}.
\]
Thus \([W^9/z^9]=c^2[e_9/z]\), since the other terms of
\((z-1)^2(z-c)^2/z\) are regular at zero. This is a generator
at the primes under consideration outside a finite set.

Apply the strip criterion with \(n=10\), \(\rho=2\), \(s=9\),
\(r=2\) and \(B=70/9\). For \(p\equiv7\bmod10\),
\(9p^2\equiv1\bmod10\), and at \(m=(9p^2-1)/10\) the
angular-constant heights are at least
\(70m/9=7p^2-7/9\) and at most \(9p^2-1\). Consequently
\[
                 9p^2\mathcal A(f^m)\equiv [z^m]P(z)^m\ne0\pmod p .
\]
There is no lower pole with even height in this strip.

For an arbitrary characteristic-zero coefficient field,
specialize the finitely generated algebra of coefficients,
preserving the prescribed face and strip, to a number field.
At all but finitely many rational primes its coefficients are
integral at primes above them. The selected scalar is defined
over the rational residue field and remains nonzero there.
Nonvanishing of each specialized moment implies nonvanishing
of the same moment before specialization. \(\square\)

There are infinitely many such fixed branch configurations,
for instance \(c\in3+341\mathbb Z\). Unordered branch sets of
distinct cross-ratios have only finitely many repeats under
Möbius transformations. Even as unmarked curves an abstract
curve of genus seven has only finitely many cyclic deck
subgroups, each determining its quotient branch set. Thus
there are infinitely many geometric curves in the construction.
No uniform density-zero exceptional set across parameters
is needed or asserted.

## Unbounded degree with a fixed primitive cover

**Corollary.** Fix an allowed branch parameter \(c\) and an integer
\(k\ge1\) prime to three. Prescribe the angular polynomial
\((z^{-1}P(z))^k\) on \(j+2i=10k\), and let all literal support obey
\[
                  j+2i\le10k,\qquad 9j+20i\ge70k.
\]
For every such fixed phase the moments of order
\(m=(9p^2-1)/(10k)\) are nonzero on a set of rational primes
\(p\equiv3^{-1}\pmod{10k}\) of natural density
\(1/\varphi(10k)\). Each phase in the strip has least sphere
polynomial degree \(11k\). The families include absolute nonpowers.

**Proof.** The progression is contained in \(p\equiv7\pmod{10}\).
Apply the strip criterion with top weight \(10k\), \(B=70k/9\),
\(s=9\), and \(r=2\). In the surviving coefficient, the exponent
of \(P\) and the selected angular index are both
\(km=(9p^2-1)/10\). Thus the computation still uses the selected
line of the same primitive order-ten factor. The progression has
the indicated density and deleting its \(\mu\)-ordinary exception
and the finite primes from coefficient specialization preserves it.

To compute the least degree, use the regular normal form with no mixed
\(x,y\) monomials. Passing to this form by \(2xy=1-t^2\) cannot increase
degree. At negative angle \(-v\), a nonzero literal coefficient
is divisible by \((1-t^2)^v\). Its height spread is therefore at least
\(2v\), whereas the strip bounds this spread by \((20k-2v)/9\). It
follows that \(v\le k\), and the greatest height \(10k+2v\) permitted
at this angle gives normal degree at most \(10k+v\).
Nonnegative angles give normal degree at most \(10k\).
The leading angular term of exponent \(-k\) has literal height
\(12k\), hence normal monomial \(y^kt^{10k}\) of degree \(11k\).

For examples, take \(F_0\) to be the regular lift of the prescribed
weight-ten polynomial, using the indicated \(-2yt^{10}\) lift for
its angle \(-1\) and monomial lifts for nonnegative angles. In the
literal chart its only weights are ten and eight. The phase
\(F_0^k+t^{10k-1}\) satisfies the strip. Its leading angular
polynomial vanishes at \(z=1\), while its weight-\((10k-1)\)
angular polynomial is the unit one. To see that it is an absolute
nonpower, work over an algebraic closure of the coefficient field
in the literal Laurent chart. Put \(z=x/t^2\), so a weight-\(w\)
layer becomes \(t^w\) times a Laurent polynomial in \(z\).
If the phase were the \(d\)-th power of a Laurent polynomial for
\(d>1\), the leading angular polynomial of that root would vanish
at \(z=1\). The next weight of its \(d\)-th power is divisible by
the \((d-1)\)-st power of that leading angular polynomial, from
the formula for the second coefficient of a power. This contradicts
the unit next layer. \(\square\)

## A complementary pair of order-five faces

The same strip method applies to the inverse deck action on
the [fixed order-five curve](CM_NEWTON_MU_ORDINARY_SELECTED_PATHS.md#a-fixed-cyclic-cover-and-a-positive-density-moment-test).
Write its equation as \(v^5=x^3(x-1)^3(x-3)^3\).
The coherent dimensions at characters \(1,2,3,4\) are now
\((1,0,2,1)\). At \(p\equiv2\bmod5\) the map from
character three to one is nonzero on the \(\mu\)-ordinary set.
The source has dimension two, so a single specified vector
does not follow from that rank statement.

Let \(A,B,C\) be an ordering of \(0,1,3\), and take the
coordinate \(u=(x-A)/(x-B)\), \(W=v(u-1)^2\).
Then \(W^5=u^3P_{A,B}(u)\), where
\[
 P_{A,B}(u)=(B-A)^6(B-C)^3 (u-1)
       \left(u-\frac{C-A}{C-B}\right)^3 .
\]
For each such cover take highest weight \(j+2i=5\)
with angular polynomial \(z^{-2}P_{A,B}(z)\), and allow all
smaller terms satisfying \(j+i\ge3\). The top regular lifts
satisfy both sides of the strip, including the negative angle
\(-2\) lift via \(y^2t^5\). The strip criterion applies with
\(s=3\), \(r=1\) and \(m=(3p-1)/5\).

To see that two fixed coordinate choices suffice at each
prime, use the original finite generator
\(e=v^3/[x(x-1)(x-3)]\) in character three. The overlap
input \([W^3/u^3]\) is represented by
\[
 e(B-A)^6 \frac{x-C}{(x-A)^2(x-B)^2} .
\]
Serre duality tests it by residues at \(A\) against
\(e^{-1}dx\) and \(xe^{-1}dx\), both regular sections of the
dualizing line times the inverse eigensheaf. For
\((A,B,C)=(0,1,3)\) the pair of residues is \((-5,-3)\);
for \((0,3,1)\) it is \((27,-81)\), with the same common
sign convention. Their determinant has absolute value
\(486=2\cdot3^5\), so they span at the relevant primes.
The nonzero directed map cannot kill both inputs. Thus at
least one of the two fixed moment tests succeeds at every
prime of the progression outside the \(\mu\)-ordinary
exception. This two-input conclusion is weaker than prescribing
one selected input in the order-ten construction.

## Exact replay

The [SageMath checker](verify_cm_newton_direct_strip_moments.py)
recomputes the primitive Frobenius polynomials at 11 and 31,
their discriminant norms, and the distinct-power checks. It
uses the power-polynomial routine in the
[order-five checker](verify_cm_newton_mu_cyclic_path.py).
For \(c=3\) it also computes the selected entries at
\(p=7,17,37\), respectively \(4,4,9\). The \(p=7\) case is
checked independently by rational expansion of the full strip
moment, including three compatible strict lower terms, and gives
normalized residue \(4\). The order-five pairing determinant
is recomputed over the rationals. These internal replays verify
the stated fixed calculations; the density and support
arguments establish the families.

    sage -python verify_cm_newton_direct_strip_moments.py
