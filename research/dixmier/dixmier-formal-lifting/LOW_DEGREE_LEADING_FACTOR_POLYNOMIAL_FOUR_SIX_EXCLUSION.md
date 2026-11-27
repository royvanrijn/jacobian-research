# Low-degree leading factors exclude polynomial orders four and six

**Theorem PolynomialFourSixLeadingFactorDegreeTwoExclusion.** Let \(k\)
have characteristic zero, and let \(P,Q\in k[x][D]\), with \([D,x]=1\),
have differential orders four and six. Suppose the leading coefficients
are \(\alpha h(x)^2\) and \(\beta h(x)^3\) with
\(\alpha,\beta\in k^\times\) and \(1\le\deg h\le2\). Then their
commutator cannot be a nonzero constant.

The leading coefficient relation is necessary for any such scalar
commutator. Its order-nine coefficient gives \(4p_4q_6'-6q_6p_4'=0\),
so \(q_6^2/p_4^3\in k^\times\). Unique factorization in \(k[x]\)
gives the displayed square and cube, up to the two constants.
The preceding [constant-leading exclusion](MONIC_POLYNOMIAL_FOUR_SIX_STRING_EXCLUSION.md)
treats \(\deg h=0\). Thus a putative polynomial order-four and
order-six string has a common leading factor of degree at least three.

## The three source substitutions

Work over an algebraic closure of \(k\). An identity over the
original field persists under this scalar extension. Scale by
nonzero constants as needed. We display the nonlinear substitutions
as embeddings in a rational Ore ring, keeping every resulting
coefficient in that ring.

If \(h\) is linear with root \(b\), set \(x=b+t^2\) and
\(D=(2t)^{-1}E\), where \(E=\partial_t\). Then \(h(x)\) is a
nonzero constant multiple of \(t^2\), and the two leading
coefficients with respect to \(E\) are constant. All new coefficients
lie in \(k[t,t^{-1}]\). The involution
\[
                  t\longmapsto -t,\qquad E\longmapsto -E
\]
fixes both original images. After scaling them to be monic, depress
the quartic by \(E\mapsto E-p_3(t)/4\), where \(p_3(t)\) is its
coefficient of \(E^3\). This coefficient is odd, so depression
commutes with the involution. In particular the new coefficient
of \(E^2\) is even in \(t\).

If \(h\) is quadratic with a double root \(b\), set \(x=b+t\),
\(E=t\partial_t\), and \(D=t^{-1}E\). The leading coefficients
are again constant and the new coefficients are Laurent polynomials.
If \(h\) has two distinct roots \(b_1,b_2\), take
\[
 b=\frac{b_1+b_2}{2},\quad d=\frac{b_2-b_1}{4},\quad
 x=b+d(t+t^{-1}),\qquad D=\frac{1}{d(t-t^{-1})}\,E,\quad E=t\partial_t .
\]
Indeed
\((x-b_1)(x-b_2)=d^2(t-t^{-1})^2\) and
\(t\,dx/dt=d(t-t^{-1})\). The leading coefficients are constant
and every new coefficient lies in
\(k[t,t^{-1},(t-t^{-1})^{-1}]\). The log derivation preserves this
localized ring. Both quadratic cases allow scaling and depression
over the same coefficient ring. Possible finite poles created by
these operations are allowed below.

We now have a monic pair satisfying
\([P,Q]=\kappa\ne0\) over \(k(t)\), with quartic subprincipal
coefficient zero and derivation
\(\delta=t^\epsilon\partial_t\), where \(\epsilon=0\) in the
linear case and \(\epsilon=1\) in both quadratic cases.

## Bounded coefficients at infinity

Use the degree at \(t=\infty\) for rational functions, numerator
degree minus denominator degree. Let \(\mathcal O_\infty\) be the
bounded functions, with maximal ideal the functions tending to zero.
Both derivations preserve this ring and ideal and induce zero
on the residue field \(k\). For \(\delta=t\partial_t\), this follows
by applying \(\delta\) to the Laurent expansion at infinity. It
annihilates the constant term and preserves all strictly negative
exponents. Finite poles play no role in this argument.

Suppose every coefficient of the already depressed \(P\) is
bounded at infinity. The monic fractional root
\(P^{1/4}=E+O(E^{-1})\) then has all coefficients in
\(\mathcal O_\infty\). The recursion for each successive coefficient
takes products and derivatives of earlier bounded coefficients
and divides by four. Thus \(B_j=(P^{j/4})_+\) belongs to
\(\mathcal O_\infty[E;\delta]\). Each \([P,B_j]\) has order at most
two, because the negative part of \(P^{j/4}\) has order at most
minus one. A residual finite operator of order \(j\) whose
commutator with \(P\) has order at most two has constant leading
coefficient, by the order-\((j+3)\) coefficient of that commutator.
This also holds at \(j=0\). Subtract a constant multiple of \(B_j\)
and descend. It follows that \(Q\) is a finite constant combination
of the \(B_j\), so is bounded as well.

Reduction now takes the pair into the commutative ring \(k[E]\)
and sends \(\kappa\) to itself, a contradiction. This tests the
coefficients after depression and needs no bound at any other
place of \(k(t)\).

## The positive highest weight

Some coefficient of the depressed \(P\) has positive degree at
infinity. For its lower coefficients \(p_i\) set
\[
 r=\max_{p_i\ne0}\ \frac{\deg_\infty p_i}{4-i}>0,\qquad
 \Delta_r=r+1-\epsilon>0 .
\]
Give \(t\) weight one and \(E\) weight \(r\), using the highest Laurent
term at infinity in each coefficient. The weighted symbol ring is
\(k[t,t^{-1},\xi]\). One Ore contraction lowers weight by at least
\(\Delta_r\), and the bracket at its first possible weight is
the Poisson bracket \(\{\xi,t\}=t^\epsilon\). The root recursion
just used also proves that \(B_j\) has weight \(jr\), since each
earlier contraction loses weight and the coefficient of the
new root term is four. Thus the order-six monic \(Q\), a constant
combination of these operators, has weight \(6r\).

The weight \(10r-\Delta_r\) is positive. For \(\epsilon=0\) this
uses \(r\ge1/4\), because a coefficient of positive integral
degree determines \(r\). For \(\epsilon=1\) it is \(9r>0\).
The leading symbols therefore Poisson-commute. The same weighted
Euler identities as in the constant-leading polynomial proof
give \(q^2=p^3\), since both symbols are monic in \(\xi\).
They have the form
\[
                    p=H^2,\qquad q=H^3,\qquad
                    H=\xi^2+\sigma t^n,\quad n=2r>0 .
\]
The symbol \(H\) is monic by factorization in the Laurent polynomial
ring. The depressed coefficient eliminates its \(\xi\) term, and
a lower coefficient contributing at the chosen maximal weight
ensures \(\sigma\ne0\). It follows that \(n\) is an integer and is
the degree of the coefficient of \(E^2\) of \(P\). In the linear
case this coefficient is even, so \(n\ge2\). In the two
quadratic cases \(n\ge1\) suffices.

Start with \(R=Q^2-P^3\), of weight less than \(12r\). The commutator
\([P,R]\) initially equals \(2\kappa Q\), of weight \(6r\).
Thus \(R\) has weight at least \(2r+\Delta_r\). If its weight
is larger, its leading symbol commutes with \(H\) for the Laurent
Poisson bracket. The weighted Euler identities force a ratio of
powers to be constant. Since \(H\) is squarefree in the Laurent
ring, the symbol is a constant multiple of \(H^j\). The possible
exponents are \(2\le j\le5\). In particular a Laurent unit \(ct^i\)
cannot multiply the power unless \(i=0\), as the ratio of powers
would then fail to be a constant.

Subtract corresponding constant multiples of \(P,Q,P^2,PQ\).
Their commutators with \(P\) have weight less than \(6r\).
Discrete weight descent therefore gives a leading symbol \(F\)
of weight \(2r+\Delta_r\) with
\[
                       \{H,F\}=\kappa H^2 .                \tag{1}
\]
For \(s^2=-\sigma\), on the parametrizations
\(t=w^2,\ \xi=\pm s w^n\) of \(H=0\), the Hamiltonian
derivation of \(H\) induces
\(\pm s w^{n-1+2\epsilon}\partial_w\). A homogeneous Laurent
polynomial of positive weight whose bracket with \(H\) is divisible
by \(H\) restricts to a constant on each parametrization, hence to
zero. The two restrictions of its degree-at-most-one remainder
in \(\xi\) force that remainder to vanish. Applying this observation
to \(F\), then to \(F/H\) of weight \(\Delta_r>0\), gives \(H^2\mid F\)
and a Laurent polynomial \(G=F/H^2\) with
\[
                            \{H,G\}=\kappa .                \tag{2}
\]

## The Laurent obstruction on a generic level

Write \(G=A(t,H)+\xi B(t,H)\) over the Laurent polynomial ring,
using \(\xi^2=H-\sigma t^n\). This representation is unique.
Let \(h_0\) stand for the indeterminate level \(H\). The polynomial
and \(\xi\) coefficients of (2) give
\[
 t^\epsilon\bigl(2(h_0-\sigma t^n)\partial_t B
               -\sigma n t^{n-1}B\bigr)=\kappa,\qquad
 \partial_t A=0 .
\]
The first identity is in \(k(h_0)[t,t^{-1}]\).
If \(B\) has a lowest exponent \(j<0\), the term
\(2h_0t^\epsilon\partial_t B\) gives a nonzero term of least exponent
\(j-1+\epsilon<0\). All contributions from the \(\sigma\) terms
have strictly larger exponent. Therefore \(B\) is a polynomial
in \(t\). It cannot be zero.

In the linear case \(\epsilon=0,\ n\ge2\).
If \(j\) is the largest exponent of \(B\), the top term on the
left has exponent \(j+n-1>0\) and nonzero coefficient
\(-\sigma(2j+n)\) times the leading coefficient of \(B\).
In the quadratic cases \(\epsilon=1,\ n\ge1\), the \(h_0\)
term contains \(t\partial_t B\), which has zero constant term.
All \(\sigma\) terms have a positive factor \(t^n\), so the
left side has zero constant term. Each case contradicts
\(\kappa\ne0\), proving the theorem.

## A local obstruction at an exact double root

**Corollary PolynomialFourSixExactDoubleRootExclusion.** Under the
same polynomial and differential-order hypotheses, a nonzero constant
commutator is impossible if the common leading factor \(h\) has a root
of multiplicity exactly two over an algebraic closure of \(k\).
The degree of \(h\) is unrestricted in this corollary.

Extend scalars to the algebraic closure and denote it again by \(k\).
At the double root \(b\), take \(\tau=x-b\).
Since \(h/\tau^2\) is a unit at \(b\), the characteristic-zero
binomial recursion gives a unit \(w\in k[[\tau]]\) with
\(h=\tau^2w(\tau)^2\). In the Laurent series field at \(t=\tau^{-1}\)
use the generator
\[
 E=-\frac{\tau w(\tau)}{w(0)}D,\qquad
 [E,t]=t\,\frac{w(t^{-1})}{w(0)} .
\]
The substitution is an injective morphism of Ore rings, since
\([D,b+t^{-1}]=1\) under the displayed change and the factor
relating \(D\) and \(E\) is nonzero. The leading coefficients
with respect to \(E\) are \(\alpha w(0)^4\) and \(\beta w(0)^6\).
Thus scaling to monic and depressing the quartic are available
over \(k((t^{-1}))\). All the resulting operators are still finite
in differential order.

The derivation on this field is
\(\delta=t(1+O(t^{-1}))\partial_t\), whose kernel is \(k\).
The bounded subring \(k[[t^{-1}]]\) and its maximal ideal are
\(\delta\)-stable, with zero residue derivation. If the depressed
\(P\) has bounded coefficients, the same fractional-root recursion
and finite leading-coefficient descent put \(Q\) in the bounded
Ore ring as well. Reduction would give a zero commutator in \(k[E]\),
contrary to the nonzero constant residue.

Otherwise choose \(r>0\) by the maximal pole degree of the lower
coefficients of \(P\), as in the quadratic log-derivation case.
The pole degrees are integers. The leading graded coefficient ring
is \(k[t,t^{-1}]\) even though the original coefficients are series.
On a term of degree \(j\) the derivation acts in that grade as
multiplication by \(j\). The \(O(t^{-1})\) term of its multiplier
only affects lower grades. Ore contractions therefore lose weight
at least \(r\) and yield the leading bracket \(\{\xi,t\}=t\)
as before. The recursion and descent for \(Q\) give weight \(6r\).
The square/cube comparison, the descent from \(Q^2-P^3\), and
equation (2) all take place in the same Laurent polynomial symbol
ring. The log-derivation obstruction on the generic level applies
verbatim and gives the contradiction.

The remaining actual polynomial case requires \(\deg h\ge3\) and
no root of multiplicity exactly two. Roots of multiplicity three or
higher, and unrestricted algebraic-function-field strings, are
outside this local obstruction.
