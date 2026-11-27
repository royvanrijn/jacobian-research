# Hecke conductor differences for the abc lifts

The two ideal-class corrections can be made over the same
abelian base. Their characters have identical inertia on
the abc support and disjoint auxiliary ramification. A
four-term conductor difference of twists of the Frey form
then detects its multiplicative primes exactly, with
the normalized weights of that base. A second abelian
base supplies the small quadratic ramification correction.
This gives an abc-equivalent height problem at a single
fixed odd level, using ordinary field discriminants and
finite conductors of twists.

## Two corrected characters over the same abelian field

Write \(R=\operatorname{rad}(abc)\), \(A=\log R\),
\(t=a/c\), and \(h(t)=\log c\) for a primitive positive
triple. Fix an odd integer \(N\ge3\) and put
\(S=\{p:p\mid R,\ p\nmid2N\}\).
Use the imaginary quadratic cover of \(S\) by distinct
fields \(J_i\). Write \(B\) for their compositum and
\(U=(B(\zeta_N))^+\). The rank bound on \(B\) gives
\([U:\mathbb Q]\le32\varphi(N)\). This extension is
abelian. If \(S\) is empty use \(\mathbb Q\) and trivial
characters throughout.
[Quadratic support cover and fixed-level correction](ABC_TOTALLY_REAL_NORM_ONE_LIFTS.md#fixed-odd-levels-with-a-single-ideal-class-correction-prime).

Fix the assigned prime ideals defining the products
\(\mathfrak a_i\) once. Apply the simultaneous class
correction twice, excluding the first repair prime
\(P_1\) from the second choice \(P_2\). Use the same
weights \(2^{i-1}\) in both products \(w_j\in B^\times\)
of norm one from \(B\) to \(B^+\).
At every place of \(B\) above \(S\) the valuations of
\(w_1,w_2\) are congruent modulo \(N\), with residue
prime to \(N\). At every other place away from \(2N\),
each \(w_j\) has valuation divisible by \(N\) except
possibly above its own correction prime. One can use a
repair prime even for a zero target class, since the
identity Artin class also admits degree-one primes.
The polynomial bound on the repair primes remains valid.

Let \(U_{\mathbb C}=B(\zeta_N)\) and take \(N\)-th roots
of \(w_1,w_2\) over \(U_{\mathbb C}\).
Complex conjugation inverts each root and \(\zeta_N\),
so it is central in the extension over \(B^+\).
Its fixed field gives a Kummer character \(\psi_j\) of
the absolute Galois group of \(U\) with values in the
complex \(N\)-th roots of unity. Use the same identification
of roots of unity for both characters. Equivalently
their restrictions to the absolute Galois group of
\(U_{\mathbb C}\) are the Kummer characters of \(w_j\).
The extensions over \(U\) are totally real, so these
finite order Hecke characters are trivial at infinity.

Put \(\psi_-=\psi_1\psi_2^{-1}\).
At every prime above \(S\), both \(\psi_j\) have inertia
of exact order \(N\) and \(\psi_-\) is unramified.
Indeed over \(U_{\mathbb C}\) this follows from the
valuation residues for \(w_j\) and \(w_1/w_2\), since
adjoining \(\zeta_N\) is unramified at \(p\nmid N\).
The quadratic extension \(U_{\mathbb C}/U\) does not
change the image of inertia under a character of odd
order. At a prime over \(P_1\), \(\psi_2\) is unramified,
so \(\psi_-\) has the same inertia as \(\psi_1\). Over
\(P_2\) it has the inverse inertia of \(\psi_2\). All
three characters are unramified at the cover prime
\(q\) when \(q\nmid 2N\), and at every remaining prime
away from \(2N\).

For example, take \(N=3\), with \(13\) assigned to
\(\mathbb Q(\sqrt{-23})\) and \(5\) to
\(\mathbb Q(\sqrt{-31})\). The repair primes \(41\)
and \(71\) allow generators
\[
 (21-2\sqrt{-23},\,9-2\sqrt{-31}),\qquad
 (30-\sqrt{-23},\,18+\sqrt{-31}),
\]
respectively. Their norms in each coordinate are the
assigned prime times the repair prime. Orient the
assigned ideals using the same square root modulo each
prime. The quotient of the two weighted norm-one products
has valuation zero at these assigned primes.
Its trace field over \(\mathbb Q(\sqrt{713})\) is
unramified at \(5,13\), and ramifies at \(41,71\).

## Exact local cancellation for the mixed twists

Let \(V\) be the Weil–Deligne system of the rational Frey
curve after base change to \(U\). One can also use its
adjoint system \(\operatorname{Ad}V\).
For a system \(T\), let \(\mathfrak f^{(2N)}(T)\)
denote its conductor ideal with the primes above \(2N\)
omitted. Define the conductor ideal quotient
\[
 \mathfrak C(T)=
 \frac{\mathfrak f^{(2N)}(T\otimes\psi_1)\,
       \mathfrak f^{(2N)}(T\otimes\psi_2^{-1})}
      {\mathfrak f^{(2N)}(T\otimes\psi_-)\,
       \mathfrak f^{(2N)}(T)}.
\]
For both systems the exact conductor identity is
\[
 \mathfrak C(V)=\mathfrak C(\operatorname{Ad}V)
       =\left(\prod_{\mathfrak p\mid p,\ p\in S}\mathfrak p\right)^2 .
                                                               \tag{TC1}
\]

At a prime over \(S\) the Frey curve is multiplicative.
The rank-two system there has unramified Weil part and
one Jordan block of size two for monodromy. More generally
\(\operatorname{Sym}^{m-1}V\) has unramified Weil part
and a single monodromy block of size \(m\), as follows
by applying the nilpotent operator to binary monomials
in characteristic zero. Thus (TC1) holds for every
symmetric power as an identity of conductor ideals.
The adjoint is the symmetric square up to an unramified
determinant twist at these places. In rank \(m\), an
unramified twist has conductor exponent \(m-1\).
A nontrivial tame inertia character kills all inertia
invariants, and the exponent is \(m\).
Thus the four exponents give \(2m-2(m-1)=2\) on \(S\).
At \(P_1,P_2\) one of the twists is unramified, and
twisting by an unramified character preserves every
conductor exponent. The terms therefore cancel in
pairs. At every other place under consideration all the
characters are unramified. These computations use the
Weil–Deligne rule \(a(r,M)=a(r)+\dim V^I-\dim(\ker M)^I\).
[Bisatt–Paterson, Lemmas 3.2–3.3](https://arxiv.org/abs/2309.02782).

The associated four-term quotient of standard
\(L\)-functions can be formed with the same two twists
in the numerator and the quotient twist and the
untwisted form in the denominator. These are twists of
one Hilbert modular form. The archimedean factors and
the field discriminant powers in their functional
equations cancel between the two sides of the quotient.
Equation (TC1) gives its remaining finite conductor
factor away from \(2N\).
For the adjoint, the conductor quotient is genuinely a
mixed-twist construction. Taking the adjoint of each
twisted rank-two form separately removes its character,
since \(\operatorname{Ad}(V\otimes\psi)=\operatorname{Ad}V\).
The required adjoint twists instead occur in the mixed
Rankin–Selberg products
\[
 L(s,(\pi\otimes\psi_i)\times
             (\widetilde\pi\otimes\psi_j^{-1}))
 =L(s,\psi_i\psi_j^{-1})\,
     L(s,\operatorname{Ad}\pi\otimes\psi_i\psi_j^{-1}),
\]
where \(\pi\) denotes the form over \(U\).
Thus comparing diagonal adjoint values of individual
twists cannot see the conductor difference.
[Artin induction and tensor formalism](https://arxiv.org/abs/1905.04282).

The same local cancellation can be expressed using three
cyclic extensions over \(U\). Let \(L_1,L_2,L_-\) be their
fields for \(\psi_1,\psi_2,\psi_-\), respectively, and let
\(\delta^{(2N)}(L/U)\) be the normalized relative log
discriminant restricted away from \(2N\). At \(p\in S\)
the two inertia indices are \(N\) and the third is one.
At either repair prime precisely the two equal inertia
costs cancel. Therefore
\[
 \frac{\delta^{(2N)}(L_1/U)+\delta^{(2N)}(L_2/U)
                       -\delta^{(2N)}(L_-/U)}{2}
   = (1-1/N)\sum_{p\in S}\frac{\log p}{e_p(U/\mathbb Q)}.
                                                               \tag{TC2}
\]
The third extension is the cyclic field attached to the
quotient character.

## A height comparison at one fixed odd level

The log norm of (TC1), divided by \(2[U:\mathbb Q]\),
equals \(\sum_{p\in S}e_p(U/\mathbb Q)^{-1}\log p\).
To recover the remaining weight, use a second
quadratic cover \(B'\) with a different exceptional
cover prime \(q'\), as in the paired construction.
If there is no exceptional prime, use \(B'=B\).
Set \(U'=(B'(\zeta_N))^+\) and let
\[
 \mathcal D(U,U')=
     d(U)+d(U')-d(UU'),\qquad d(F)=\log\operatorname{rd}(F).
\]
Away from \(2N\), taking the maximal real subfield
preserves the inertia of \(B(\zeta_N)\): inertia acts
trivially on \(\zeta_N\) while complex conjugation
does not. At every \(p\in S\), the two abelian bases
therefore have identical inertia indices \(e_p\in\{1,2\}\).
Their differing cover primes ramify in at most one.
Tame inertia in the compositum has order the least
common multiple, and the omitted normalized different
is \(O(\log(2N))\) by (PD2) and the degree bounds.
It follows that
\[
 \mathcal D(U,U')+
  \frac{\log\mathrm N_{U/\mathbb Q}\mathfrak C(V)}
       {2[U:\mathbb Q]}
     = A+O(\log(2N)).
                                                               \tag{TC3}
\]
One may replace \(V\) by \(\operatorname{Ad}V\).
If all three field discriminants are also restricted
away from \(2N\), the left side is exactly
\(\sum_{p\in S}\log p\). The formula retains the whole
tame radical even when \(N=3\) is fixed once and for all.
[Paired cover choices and normalized different](ABC_TOTALLY_REAL_NORM_ONE_LIFTS.md#cancelling-the-auxiliary-discriminants-by-a-second-lift).

Full abc is equivalent to an upper bound for \(h(t)\)
by \(1+\varepsilon\) times the left side of (TC3), plus a
constant depending on \(\varepsilon\) and this one fixed
\(N\), uniform in the triples. All fields in the base
discriminant terms are abelian, and the four twists are
over the same field. The open step is a comparison of
integral periods or congruence invariants for the mixed
twists which keeps the four-term conductor difference.
Tilouine and Urban relate integral periods and
congruence ideals for quadratic base change at a
residue prime under hypotheses on the localized
cohomology. Their power-of-two degree period comparison
applies to balanced periods under additional minimality
hypotheses. At \(N=3\), the abelian bases here do have
power-of-two degree. The comparison concerns the
untwisted base change, and an estimate for the mixed
odd-order twist periods is still required.
Pasten's rational height estimate via congruence
multiplicities gives a bound in terms of
\(M\varphi(D)\log\ell\) at an admissible rational
factorization \(DM\) of the conductor. The totally
real estimate in Theorem 18.10 also has a
parametrization-degree term. Neither estimate supplies
the required relative comparison.
[Tilouine–Urban, Theorems 4.5 and 4.28](https://arxiv.org/abs/1811.11166).
[Pasten, Theorem 18.10](https://arxiv.org/abs/1705.09251).
[Pasten, Theorem 1.1](https://arxiv.org/abs/2609.17390).

A fixed cubic Dirichlet twist over the rationals also gives a
single adjoint value whose completed size recovers the radical.
The associated normalization by the Néron period gives a sharp
height comparison for algebraic values in one fixed quadratic field.
[Fixed cubic adjoint value over the rationals](ABC_FIXED_CUBIC_ADJOINT_VALUE.md).
