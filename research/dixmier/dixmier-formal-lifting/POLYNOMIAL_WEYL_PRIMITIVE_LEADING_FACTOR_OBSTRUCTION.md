# A high-root obstruction for primitive Weyl leading factors

Let \(P,Q\) be complex polynomial differential operators of orders
\(ag,bg\), with \(a,b\ge2\) coprime, and let \([Q,P]=\kappa\ne0\).
Write their leading coefficients as scalar multiples of \(h^a,h^b\).
For nonconstant \(h\) let \(e\) be the gcd of \(g\) and all the
geometric root multiplicities of \(h\).

## A resultant with a controlled differential leading coefficient

**Proposition PolynomialWeylPrimitiveLeadingResultant.** Suppose
\(g\ge2\) and \(e=1\), and normalize the leading coefficient of
\(P\) to \(h^a\). Write \(U=P^{1/a}\) for the rational
pseudodifferential root with leading coefficient \(h\).
There is a Laurent polynomial \(G\in\mathbb C[U,U^{-1}]\)
of degree \(b\) with leading coefficient \(\gamma\ne0\) such
that
\[
              T=G(U)-Q,\qquad
              \operatorname{ord}T=1-ag,\qquad
              \operatorname{lc}T=\frac{\kappa B}{ag\,h^a},
              \quad B'-\frac{h'}{gh}B=1 .
\]
For some integer \(M\ge0\), the commutative polynomial
\[
  \mathcal H(p,q)=p^M\prod_{\omega^a=1}(q-G(\omega u)),\qquad u^a=p
\]
has top \((a,b)\)-weight expression \(p^M(q^a-\gamma^a p^b)\).
With coefficients in \(P\) placed on the left, the ordinary
operator \(W=\mathcal H(P,Q)\) satisfies
\[
 \begin{split}
  D &= aM+(a-1)b-a,\\
  \operatorname{ord}W &= Dg+1,\qquad
  \operatorname{lc}W=-\frac{\kappa\gamma^{a-1}}{g}B h^D,\qquad
  [W,P]=\kappa\mathcal H_q(P,Q).
 \end{split}                                                   \tag{1}
\]

The fixed-point part of the [spectral phase calculation](POLYNOMIAL_WEYL_SMALL_ROOT_SPECTRAL_CONSTRAINTS.md#common-spectral-phases-and-the-opposite-fourier-coefficient)
does not require a small root. Indeed on the cover containing
\(p_n^{1/(ag)}\), where \(p_n=h^a\), spatial diagonalization has
simple leading eigenvalues \(\omega_i p_n^{-1/(ag)}Z\),
\(Z^{ag}=\lambda\). Flatness makes the phase coefficients above
degree one constant in \(x\). The spatial deck on the leading
root has order \(g/e=g\), so all those phase degrees are
multiples of \(g\). The degree-one term contributes only at
negative differential order in \(Q\). The [formal wave
calculation](POLYNOMIAL_WEYL_SMALL_ROOT_SPECTRAL_CONSTRAINTS.md#fractional-powers-in-the-differential-mate)
then gives \(Q=(F(U))_+\) for a polynomial \(F\) of degree \(b\)
with constant coefficients. Its leading coefficient is denoted
by \(\gamma\). This uses the growing spectral degree at an
ordinary point, even if every zero of \(h\) is high.

Start with \(F(U)-Q\), which has negative order and
commutator \(-\kappa\) with \(P\). At any order \(r>1-ag\)
of this operator or of a remainder obtained by subtracting
commuting terms, the top commutator coefficient forces
its leading coefficient to be proportional to \(h^{r/g}\).
Rationality and \(e=1\) require \(g\mid r\).
Subtract that multiple of \(U^{r/g}\) and continue until
order \(1-ag\). The nonzero scalar commutator prevents
the order from dropping below \(1-ag\).
The coefficient at that order solves
\((1-ag)t p_n'-ag p_n t'=-\kappa\).
The displayed formula for \(\operatorname{lc}T\) follows
from the [polynomial residue primitive](POLYNOMIAL_WEYL_LEADING_FACTOR_RESIDUE.md).
Its difference from any other solution is proportional to
\(h^{(1-ag)/g}\), which is not rational when \(e=1\) and \(g\ge2\).

The product defining \(\mathcal H\) is invariant under
\(u\mapsto\omega u\), so its coefficients lie in
\(\mathbb C[p,p^{-1}]\) before multiplication by \(p^M\).
Coprimality of \(a,b\) makes the top product
\(q^a-\gamma^a p^b\). In evaluating \(W\), substitute
\(Q=G(U)-T\) into the ordered powers of \(Q\).
The terms with no \(T\) sum to zero. Each \(T\) lowers
differential order by \((a+b)g-1>0\) from that of \(G(U)\).
In the terms linear in \(T\), moving it past a coefficient
or a \(G(U)\) factor costs at least one additional
differential order. Thus the leading coefficient is
that of the commuting linear term
\(-\mathcal H_q(P,G(U))T\). Its leading coefficient is nonzero
because the top \(q\)-derivative is \(a p^M q^{a-1}\).
This gives the order and leading coefficient in (1).
The last identity in (1) is the ordinary product rule
using \([Q,P]=\kappa\).

## No high root in the primitive case

**Theorem PolynomialWeylPrimitiveLeadingFactorSingleSmallRoot.**
Under the above hypotheses with \(e=1\) and \(g\ge2\), the common
leading factor \(h\) has a single geometric root, of multiplicity
strictly smaller than \(g\). In particular every multi-root common
leading factor in a nondividing pair of differential orders with
\(g\ge2\) has a multiplicity divisor \(e>1\) common to all roots
and to \(g\).

Suppose a root of multiplicity \(\mu>g\) exists and translate
it to zero. Start at differential order and then increase the
position penalty of weight \((-s,1)\) from zero.
For small positive \(s\), the maximizing points are
\(a(\mu,g),b(\mu,g)\).
If neither pair has a first edge before \(s=g/\mu<1\),
both weights at \(s=g/\mu\) are zero, contrary to the
scalar commutator and positive contraction loss.
At the first edge the same two points are the endpoints
of largest differential exponent. In Fourier coordinates
these are the aligned endpoints in
[Guccione–Guccione–Valqui, Lemma 7.1](sources/guccione-guccione-valqui-1111.6100v3.pdf).
The lemma gives positive weights, zero leading bracket, and
simultaneous edges. Thus the top forms are \(R^a,\gamma R^b\)
with \(R\) polynomial and of differential degree \(g\):
the positive weight ratio and the zero bracket give equal
powers by the homogeneous Euler identities, and coprime
factorization gives \(R\).

Rescale the edge weight to \((-1,\rho)\), \(\rho>1\).
Put \(w=\rho g-\mu>0,\ L=\rho-1>0\) and \(z=x^\rho y\);
a common ramified \(x\)-cover can be used in this notation.
Then
\[
                            R=x^{-w}f(z),\qquad \deg f=g .
\]
The constant term of \(f\) vanishes because \(R\) is
polynomial and \(w>0\). The edge has more than one
term, so \(f\) has at least two distinct roots, including zero.

The top monomial given by (1) gives \(W\) weight at least
\(Dw+L\), since \(B\) has a simple zero at the translated
root. The top expression of \(\mathcal H_q(P,Q)\) has
weight \((D+a)w\) and is \(a\gamma^{a-1}R^{D+a}\).
If \(W\) had weight \(v>Dw+L\), its top bracket with \(R^a\)
would vanish. Writing its top as \(x^{-v}A(z)\), one
would have \(A\) proportional to \(f^{v/w}\).
If \(d_f\) is the gcd of the root multiplicities of \(f\),
then \(d_f v/w\) is an integer. Since \(f\) has at least
two distinct roots, \(g/d_f\ge2\). Thus
\(\deg A=gv/w\ge Dg+2\), contrary to (1).

The weight of \(W\) is consequently \(Dw+L\).
Write its top as \(x^{-(Dw+L)}A(z)\), \(\deg A\le Dg+1\).
The top equation from (1), after dividing by \(R^{a-1}\),
is \(\{\operatorname{in}W,R\}=\kappa\gamma^{a-1}R^{D+1}\).
In the variable \(z\) it reads
\[
            (Dw+L)A f'-w A'f=\kappa\gamma^{a-1} f^{D+1}.
\]
At a root of \(f\) of multiplicity \(\nu\), if
\(\operatorname{ord} A\le D\nu\) then the left leading
coefficient cannot vanish because \(L>0\). Its order
would be at most \((D+1)\nu-1\), below that on the right.
Thus \(A/f^D\) is a polynomial of degree at most one
and vanishes at every root of \(f\). There are at
least two roots, a contradiction.

The [leading-factor residue theorem](POLYNOMIAL_WEYL_LEADING_FACTOR_RESIDUE.md)
forces a high root if there are multiple roots and
excludes multiplicity exactly \(g\). This proves the result.

## Consequences for centered rectangles and low orders

**Corollary PolynomialWeylCenteredCornerNoncoprime.** In a centered
one-pole rectangle with normalized corner \((u,v)\), \(0<u<v\),
and nondividing ratios \(a:b\), one has \(\gcd(u,v)>1\).
More precisely the opposite Fourier leading factor has gcd
of its multiplicities with \(u\) at least two.

In the opposite Fourier chart the orders are \(au,bu\)
and the common leading factor has degree \(v\), from the
two original corner monomials. The reduced-order bound
already excludes \(u=1\). If the multiplicity divisor
in that chart were one, the preceding theorem would give
a single root of multiplicity below \(u\), contrary
to \(v>u\). The multiplicity divisor also divides \(v\).

**Corollary PolynomialWeylOrderTenAndElevenGeneration.** Over
characteristic zero, a polynomial Weyl pair with nonzero scalar
commutator generates if one member has differential order ten
or eleven.

For order ten the [fixed-member reduction](POLYNOMIAL_WEYL_TRANSLATED_ROOT_CROSSINGS.md#the-order-ten-case)
leaves orders \(10,15\), \(g=5\), with a high root of
multiplicity twenty and small roots of multiplicity below
five. Since five is prime, the multiplicity divisor is one,
contradicting the theorem. For order eleven subtract powers
of the fixed member while the other order is divisible by
eleven. The top coefficient comparison permits each subtraction
and the order drops. The other order cannot be a proper divisor
of eleven in a proper pair by the reduced-order bounds.
The residual gcd would be one, already excluded.
For a general characteristic-zero field the finite coefficient
field of a pair embeds in \(\mathbb C\). Generation descends
since expressions for the generators solve finite linear
systems over that coefficient field.

**Corollary PolynomialWeylProperOrderTwelveBound.** In a proper
characteristic-zero polynomial Weyl pair with nonzero scalar commutator,
each differential order is composite and at least twelve.

The [order-eight generation theorem](POLYNOMIAL_WEYL_TRANSLATED_ROOT_CROSSINGS.md)
already gives composite orders at least ten. The preceding corollary
excludes ten and eleven, giving the stated lower bound. The earlier
lower bound remains valid independently of the new primitive-factor
argument.

**Corollary PolynomialWeylOrderTwelveCommonMultiplicityDivisor.**
In the existing differential-order-twelve residual list,
the multiplicity divisor of the nonconstant common leading
factor is greater than one. In particular the four-simple-root
degree-thirty \((12,18)\) profile \(A\,x^{26}(x^4-\lambda_0)\)
is excluded. In the branch with high multiplicity 26 at
degree thirty, the remaining multiplicities can only be
\((4)\) or \((2,2)\).

Each row of the [order-twelve residual list](POLYNOMIAL_WEYL_ORDER_TWELVE_RESIDUAL.md)
has a high root, so the primitive case is impossible.
With multiplicity 26 and gcd six the divisor must be two,
and the remaining degree four is partitioned into even
multiplicities. The [even-profile analysis](POLYNOMIAL_WEYL_ORDER_TWELVE_HIGH_ROOT_FRACTIONAL_EXCLUSION.md#the-remaining-even-profiles)
imposes further local restrictions on those partitions.
It also gives an actual order-nine polynomial operator \(W\) whose
leading factor with \(P\) is primitive for the new gcd three. The
[length-two commutator reduction](POLYNOMIAL_WEYL_ORDER_TWELVE_HIGH_ROOT_FRACTIONAL_EXCLUSION.md#a-primitive-length-two-commutator-problem)
retains the defining identity for \(W\) and computes its first
order-34 binomial resultant. The new pair has a nonzero critical-face
bracket, so the scalar theorem above does not apply to it. Excluding
the two resulting even profiles remains open.

## Exact controls

The [checker](verify_weyl_primitive_leading_resultant.py) uses
only the Python standard library. From the repository root run
`python3 research/dixmier/dixmier-formal-lifting/verify_weyl_primitive_leading_resultant.py`.
It computes 50 commutative Laurent resultants by a multiplication
matrix, checks their top weighted terms and the leading derivative
after substitution, and checks three exact multiplicity-26 face
identities and the remaining degree-thirty partitions. It also
checks the universal length-two chain identities, a corrupted
defining-identity control, and the nonzero order-34 resultant faces.
These are finite arithmetic controls for the written proofs; they
do not exclude the even profiles.
