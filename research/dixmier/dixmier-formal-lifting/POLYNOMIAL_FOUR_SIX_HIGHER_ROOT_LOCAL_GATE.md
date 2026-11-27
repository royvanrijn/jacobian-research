# Higher-root constraints for polynomial orders four and six

**Theorem PolynomialFourSixHigherRootLocalGate.** Let \(k\) have
characteristic zero, let \([D,x]=1\), and let \(P,Q\in k[x][D]\)
have orders four and six with \([P,Q]=\kappa\in k^\times\).
Write their leading coefficients as \(\alpha h^2,\beta h^3\)
with \(\alpha,\beta\in k^\times\). At any root \(b\) of \(h\) over
\(\overline k\) of multiplicity \(m\ge3\), one has \(m\equiv2
\pmod4\). In the local monic normalization with derivative proportional
to \(\sqrt h\,D\), depress the quartic by shifting its cubic coefficient
to zero. Let \(r\) be the largest pole order of a lower coefficient
divided by its deficit in differential order. Then \(r\) is a positive
integer and exactly one of the following necessary alternatives holds:
\[
          m=4r+2,\qquad\text{or}\qquad m=8r+2\ \text{with }r\ge2 .
\]
The original coefficients \(p_i,q_i\) of \(D^i\) are divisible by
\((x-b)^i\) for \(i\ge2\), and
\[
                  p_1(b)q_0'(b)-q_1(b)p_0'(b)=\kappa.
\]

The leading-coefficient form follows from the order-nine commutator
identity \(4p_4q_6'-6q_6p_4'=0\). The
[exact double-root exclusion](LOW_DEGREE_LEADING_FACTOR_POLYNOMIAL_FOUR_SIX_EXCLUSION.md#a-local-obstruction-at-an-exact-double-root)
eliminates multiplicity two by a separate case. Simple roots are not
restricted here.

## Local normalization and the leading symbol

Work over the algebraic closure of the coefficient field and continue
to write \(k\). Put \(\tau=x-b=t^{-e}\), where \(e=1\) for even \(m\)
and \(e=2\) for odd \(m\). A square root \(\lambda\) of \(h\) lies
in \(k((t^{-1}))\), since \(em\) is even. Take \(E\) to be a nonzero
constant multiple of \(\lambda D\). This is a change of derivative
generator in a Laurent series Ore ring, with
\[
    [E,t]=t^{1-s}(1+O(t^{-1})),\qquad s=e(m-2)/2>0,
\]
after choosing the constant multiple to make the leading multiplier one.
Scale the two operators to be monic and denote their nonzero scalar
commutator by \(\kappa_*\). Their degree-zero coefficients in \(E\)
remain regular at \(\tau=0\): every positive power of \(\lambda D\)
or of its inverse-scaled derivative has a right factor equal to the
corresponding derivative, so contributes nothing to the constant
coefficient. Write the monic quartic and sextic before depression as
\(P_0,Q_0\). Set \(a\) equal to one fourth the cubic coefficient of
\(P_0\) and use \(E\mapsto E-a\) to get the depressed pair
\(\widetilde P,\widetilde Q\). This shear is an Ore automorphism.
The reverse shear recovers \(P_0,Q_0\).

For the derivation induced by \(E\), the pole degree in \(t\) of a
derivative goes down by at least \(s\). In particular the bounded
coefficient ring and its maximal ideal are stable, with zero residue
derivation. If the coefficients of \(\widetilde P\) were bounded, the
same would hold for \(\widetilde Q\). Indeed the monic formal fourth
root \(L\) of \(\widetilde P\) has bounded coefficients by the
coefficient recursion for \(L^4\). Each finite part \(B_j=(L^j)_+\)
has commutator of order at most two with \(\widetilde P\), since the
negative part has order at most minus one. An order-\(j\) operator with
commutator of order at most two with \(\widetilde P\) has constant
leading coefficient, from the order-\((j+3)\) coefficient. The
derivation has kernel \(k\). Descent by the \(B_j\), including \(B_0=1\),
expresses \(\widetilde Q\) as a constant combination of them.
Reduction in the bounded ring would then give a zero commutator in
\(k[E]\), contrary to \(\kappa_*\).

Thus \(r>0\). Assign weights \(\operatorname{wt}(t)=1\) and
\(\operatorname{wt}(E)=r\). In the associated Laurent polynomial
symbol ring \(k[t,t^{-1},\xi]\) each product contraction loses at
least \(r+s\), and the first possible commutator symbol uses
\(\{\xi,t\}=t^{1-s}\). The root recursion gives \(B_j\) weight \(jr\):
comparison in \(L^4\) at each new root coefficient bounds its term
by weight \(r\), since contractions have positive loss. The term
\(E^j\) realizes \(jr\). Thus the two depressed operators have
weights \(4r,6r\). The scalar commutator forces \(9r\ge s\).

In the original monic pair the reverse shear recovers the bounded
constant coefficients. If \(\deg_t a=d>r\), weight \(d\) on the
derivative gives the reverse-sheared quartic constant coefficient
leading term \(a^4\): every contraction loses positive weight and
the terms coming from the lower coefficients have weight
at most \((4-i)r+id<4d\). This contradicts boundedness.
Hence \(\deg_t a\le r\). Over a fractional Laurent extension write
the two leading symbols as \(t^{4r}f(\eta), t^{6r}g(\eta)\),
where \(\eta=\xi/t^r\) and \(f,g\in k[\eta]\) are monic.
The reverse shear translates \(\eta\) by a constant \(c\), which
is the leading coefficient of \(a\) when \(\deg_t a=r\) and zero
otherwise. Boundedness of the original constant coefficients therefore
gives \(f(c)=g(c)=0\).

If \(9r=s\), the leading commutator identity is
\[
                    2r(3f'g-2fg')=\kappa_* .
\]
Evaluation at \(c\) contradicts the common root. In particular
\(9r>s\) and the two leading symbols Poisson-commute. The same
identity with right side zero gives \(g^2=f^3\) by monicity.
The depressed cubic coefficient gives
\[
    p=H^2,\qquad q=H^3,\qquad H=\xi^2+\sigma t^n,
                \qquad n=2r,\quad \sigma\ne0 .
\]
Here \(p,q\) denote the leading symbols of \(\widetilde P,\widetilde Q\).
The nonzero lower coefficient that attains \(r\) ensures
\(\sigma\ne0\), and the \(\xi^2\) coefficient makes \(n\) integral.
The common root \(c\) satisfies \(c^2=-\sigma\), so \(c\ne0\).
Thus \(r=\deg_t a\) is an integer. After reversing the shear
the common quadratic has the new leading form
\(\xi(\xi+2c t^r)\).

## The Hamiltonian equation on a generic level

Start with \(R=\widetilde Q^2-\widetilde P^3\).
Its weight is less than \(12r\), and
\([\widetilde P,R]=2\kappa_*\widetilde Q\) has weight \(6r\).
As long as this leading commutator is preserved, the weight of
\(R\) is at least \(3r+s\). At a strictly larger weight its
leading symbol commutes with \(H\). For positive weighted
homogeneous symbols, the Euler identities and this bracket
force a ratio of powers with \(H\) to be constant. Since
\(H\) is squarefree, such a symbol is a constant multiple
of an integral power \(H^i\). Its weight bounds require
\(2\le i\le5\). Subtract the corresponding constant multiple
of \(\widetilde P,\widetilde Q,\widetilde P^2,\widetilde P\widetilde Q\).
Their commutators have lower weight than \(6r\), and descent
terminates at weight \(3r+s\). The resulting leading symbol \(F\)
satisfies \(\{H,F\}=\kappa_* H^2\).

On the parametrizations \(t=w^2,\ \xi=\pm\sqrt{-\sigma}\,w^n\) of
\(H=0\), the Hamiltonian derivation induces a nonzero monomial multiple
of \(d/dw\). Any Laurent homogeneous symbol of positive weight
whose bracket is divisible by \(H\) restricts to zero on both
parametrizations, since its restrictions must be constants of
positive degree. Division by the monic quadratic in \(\xi\)
then shows divisibility by \(H\). Apply this first to \(F\)
and then to \(F/H\), which still has positive weight. It
follows that \(G=F/H^2\) is a Laurent polynomial and
\(\{H,G\}=\kappa_*\).

Write \(G=A(t,H)+\xi B(t,H)\) and view \(T=H\) as the indeterminate
level. The equation for the part of \(\xi\)-degree zero is
\[
      2T B_t-\sigma t^{n-1}(2t B_t+nB)=\kappa_* t^{s-1}
                        \quad\text{in }k(T)[t,t^{-1}].
\]
A negative least exponent of \(B\) is impossible: its derivative
contributes a nonzero term of strictly smaller exponent than
all terms from the \(\sigma\) part and the right side. For each
residue class of exponents modulo \(n\), the largest exponent
contributes a nonzero highest term through the \(\sigma\) part,
since the exponent is nonnegative. Consequently only the
residue class of \(s\) can occur. The largest exponent in
\(B\) is \(s-n\), and its least exponent must be zero. Otherwise
its derivative gives a nonzero term of exponent less than
\(s-1\) that cannot cancel. Therefore \(s=nJ\) with positive
integer \(J\), and \(9r>s\) gives \(J\le4\).
Since \(n=2r\) is even, \(s\) is even. This rules out odd \(m\),
where \(s=m-2\), and in the even case gives \(m=4rJ+2\).

## The original coefficients and the maximal ideal

Now \(m\) is even, so \(e=1\). By induction on \(j\) the coefficient
of \(D^i\) in \(E^j\) has \(\tau\)-valuation at least
\(i+j(m/2-1)\). Applying \(E\) on the left either differentiates
the coefficient and multiplies by \(\lambda\), or gives one
additional \(D\) and multiplies by \(\lambda\).
The reverse-sheared monic operators retain bounds \(4r,6r\) for
their weights. Translating from their \(E\)-coefficients to the
original \(D\)-coefficients, the smallest lower bound occurs
at \(j=i\), since \(m/2-1+r>0\). It follows that
\[
 \operatorname{ord}_\tau p_i\ge mi/2-(4-i)r,\qquad
 \operatorname{ord}_\tau q_i\ge mi/2-(6-i)r .
\]
Both inequalities are strict at \(i=1\), from the leading
symbols \(\xi^2(\xi+2ct^r)^2,\xi^3(\xi+2ct^r)^3\)
of the undepressed pair. At \(i\ge2\) both lower bounds
are at least \(i\), using \(r\le(m-2)/4\).

If \(J\ge3\), or \(J=2\) with \(r=1\), then \(m/2\ge5r\).
The strict estimates at \(i=1\) make both first-order
coefficients vanish at \(b\). All terms of both operators
preserve \((\tau)\subset k[[\tau]]\), so the commutator acts
by zero on its one-dimensional quotient. This contradicts
\(\kappa\ne0\), leaving the two claimed alternatives.
Finally apply \([P,Q]\) to \(1\) and evaluate at \(b\).
All terms of order at least two vanish and the degree-zero
products cancel, giving the asserted first-order identity.
