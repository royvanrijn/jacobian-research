# Exclusion of constant-leading Weyl counterexamples at order thirty-two

**Theorem MonicPolynomialOrderThirtyTwoGeneration.** Let \(k\) have
characteristic zero and let \(P,Q\in A_1(k)\) have nonzero scalar
commutator. If one member has differential order \(32\) and nonzero
constant leading coefficient, the two operators generate \(A_1(k)\).

Together with the
[order-thirty-one theorem](POLYNOMIAL_WEYL_WEIGHTED_DESCENT_AND_ORDER_THIRTY_ONE.md),
this gives generation for every constant-leading member of order at
most \(33\). The argument at order \(32\) uses the four faces
from [MonicPolynomialOrderThirtyTwoFaces](MONIC_POLYNOMIAL_ORDER_THIRTY_TWO_FACES.md).
The root shifts below take place in explicitly specified fractional
Weyl algebras.

## Successors after a root shift

Suppose there is a proper pair, and extend constants to an algebraic
closure. Use the Fourier presentation \([Y,X]=1\) of the four-face
theorem and scale the mate to \([Q,P]=1\). The face direction is
\((\ell,-1)\), with \(\ell=4\) or \(2\), and its normalized common
factor is
\[
                 R=X L(XY^\ell).
\]
The highest \(Y\)-exponent corner of \(R\) on this face is
\((4,12)\) or \((6,10)\). In the variable \(\eta=X^{1/\ell}Y\),
the polynomial \(L(\eta^\ell)\) has maximal root multiplicity
\(\mu\), where the possibilities are
\[
                (\ell,\mu)=(4,2),(4,3),(2,4),(2,5).
\]
All these roots are nonzero. The face of \(P\) is proportional
to \(R^2\) and that of \(Q\) to \(R^b\), with \(b\ge3\) odd.

Use [Guccione–Guccione–Valqui II, Proposition 4.3, printed pages 24–26](https://arxiv.org/pdf/1205.6827v1).
Their fractional Weyl ring \(W(l)\) is \(k[X^{1/l},X^{-1/l}][Y]\)
with the derivation \(d/dX\). A shift at a direction
\((\rho,\sigma)\) adjoins fractional powers with denominator
\(\operatorname{lcm}(l,\rho)\) and sends
\(Y\) to \(Y+\lambda X^{\sigma/\rho}\). Here \(\lambda\) is a root
of maximal multiplicity in the polynomial on the leading face
of \(P\), including a possible factor from its least \(Y\)-exponent.
The shift preserves the scalar commutator.

Initially the hypotheses of Proposition 4.3 hold as follows.
The second coordinate of the direction is negative, and
\([Q,P]=1\). Both faces have multiple terms, their positive
weights have ratio \(2:b\), and their common powers give zero
leading bracket. The \((1,-1)\) weights of their two highest
\(Y\)-exponent endpoints are negative because the Bernstein
corners have first coordinate smaller than the second.
Thus the proposition applies to a root of multiplicity \(2\mu\)
in the face of \(P\). After the first shift, its successor
direction is \((\rho',\sigma')\), written as weights \((1,-s)\)
by positive scaling, and both operators are in \(W(\ell)\).
The successor has
\[
 C=\left(1+\frac{\mu}{\ell},\mu\right),\qquad
 C_P=2C,\qquad C_Q=bC,\qquad
             \frac1\ell<s<\frac{C_x}{\mu}<1,              \tag{1}
\]
where \(C_P,C_Q\) are the endpoints of largest \(Y\)-exponent.
Proposition 4.3(6)–(7) gives these endpoints after the
shift. Items (1), (3), (4) give the strict increase of the
slope \(s\), positive weights \(2w,bw\) with \(w=C_x-\mu s\),
and multiple terms in both faces.

At every step below that has zero leading bracket we use
Proposition 4.3 again. In such a step the scalar commutator
has been preserved, the direction still has negative second
coordinate and both faces have multiple terms and positive
weights in ratio \(2:b\). Their highest \(Y\)-exponent endpoints
have negative \((1,-1)\) weights, as checked for each corner
below. These verify conditions (a)–(g) of the proposition.
It provides another edge of strictly larger slope, positive
weights of the same ratio, and the endpoint formulas (6)–(7).

## The bracket at the first successor cannot be scalar

At a successor with endpoint \(C\) in (1) the highest possible
commutator weight is \((b+2)w-(1-s)\ge0\). If the leading
bracket is nonzero, this weight must be zero. Put \(B=b+2\);
then
\[
   \ell s=\frac{B(\ell+\mu)-\ell}{B\mu-1}.
\]
Set \(\beta=\mu(\ell-1)-\ell>0\). The gcd of numerator and
denominator in this fraction divides \(\beta\), as is seen by
subtracting \((\ell+\mu)\) times the denominator from \(\mu\)
times the numerator. Since the \(X\)-exponents are in
\(\ell^{-1}\mathbb Z\), every deficit from the top \(Y\)-exponent
\(2\mu\) in the face of \(P\) is divisible by the reduced
denominator of \(\ell s\).

In the cases \((\ell,\mu)=(4,2),(2,4)\), \(\beta=2\) and
\(B\mu-1\) is odd, so the denominator is larger than \(2\mu\).
For \((4,3)\), \(\beta=5\). The only possible denominator at
most \(6\) occurs at \(B=7\) and equals \(4\).
Indeed \(5\mid(3B-1)\) requires \(B\equiv2\bmod5\), and the
next odd \(B\) with that property is \(17\).
For \((2,5)\), \(\beta=3\). The only possible denominator
at most \(10\) occurs at \(B=5\) and equals \(8\).
When the denominator exceeds \(2\mu\) the first leading face
is a monomial, contrary to edge membership. In the two exceptional
cases its polynomial in \(\eta=YX^s\) is divisible by
\(\eta^2\), since \(2\mu\) is congruent to \(2\) modulo the
denominator. Its Poisson bracket with the mate then vanishes
at \(\eta=0\) and cannot equal a nonzero scalar.
This excludes a scalar leading bracket at every successor
with the same \(\ell,C\) and slopes in (1).

## The coefficient equation and the lattice table

At a commuting successor with endpoint \(C\), the leading
symbols are proportional to \(R_1^2,R_1^b\), where
\[
                      R_1=X^w f(\eta),\qquad \eta=YX^s,
                      \qquad \deg f=\mu.
\]
This is the coprime power comparison
[Guccione–Guccione–Valqui II, Theorem 2.12(2)](https://arxiv.org/pdf/1205.6827v1).
In particular \(R_1\) belongs to the Laurent symbol ring
\(k[X^{1/\ell},X^{-1/\ell},Y]\).
If \(f\) has a single distinct root, it is nonzero since the
face is not a monomial. Adjacent \(Y\)-exponents occur in
\(R_1\), forcing \(\ell s\) integral. The root shift of
Proposition 4.3 then keeps the same denominator \(\ell\)
and the same normalized highest \(Y\)-endpoint \(C\), because
the chosen root has multiplicity \(\mu\).
The intervals (1) allow no such slope for \(\ell=2\), and
only \(s=1/2\) for \(\ell=4\). Such a shift can therefore
occur at most once while the corner remains \(C\).
In all four cases \(C_x-\mu<0\), as required for the
new application of Proposition 4.3.

Now suppose \(f\) has at least two distinct roots. The commuting
edge satisfies all of Proposition 4.3(a)–(g), as checked above.
Its item (11) gives a homogeneous \(F\) of weight \(1-s\)
with \([P,F]_{\rho',\sigma'}=\ell_{\rho',\sigma'}(P)\).
Dividing the resulting Poisson equation by a scalar multiple
of \(2R_1\) gives \(\{R_1,F\}=cR_1\), with \(c\ne0\).
Write the symbol as \(F=X^{1-s}W(\eta)\). The bracket is
\(X^w((1-s)f'W-w fW')\), up to sign. Therefore
\[
                 W'-\frac{1-s}{w}\frac{f'}{f}W=c_0,
                 \qquad c_0\ne0.                         \tag{2}
\]
A finite pole of \(W\) is impossible because its leading
negative exponent minus a nonnegative multiple of the
root multiplicity is nonzero. Thus \(W\) is polynomial,
with a simple zero at each root of \(f\).
A root of \(f\) with multiplicity \(m_0\) is impossible
if \((1-s)m_0/w=1\). Since \(f\) has at least two roots,
comparison at infinity gives
\[
 M=\deg W=\frac{\mu(1-s)}w\ge2,\qquad
             \ell s=\ell-\frac{\beta M}{\mu(M-1)}.       \tag{3}
\]
Let \(j<\mu\) be the least power of \(\eta\) in \(f\).
The reduced denominator \(h\) of \(\ell s\) divides
\(\mu-j\), and \(f=\eta^j H(\eta^h)\) with \(H(0)\ne0\).
Also \(M-1\mid\beta(\mu-j)\) by (3) and
\(\gcd(M,M-1)=1\). Therefore
\(2\le M\le\beta\mu+1\), giving the following exact table
from (1) and (3). The listed \(j\)'s are the integers
in \([0,\mu-1]\) with \(h\mid\mu-j\).

| \((\ell,\mu)\) | \(M\) | \(s\) | \(h\) | \(j\) |
|---|---:|---:|---:|---|
| \((4,2)\) | 2 | \(1/2\) | 1 | \(0,1\) |
| \((4,2)\) | 3 | \(5/8\) | 2 | \(0\) |
| \((4,3)\) | 3 | \(3/8\) | 2 | \(1\) |
| \((4,3)\) | 6 | \(1/2\) | 1 | \(0,1,2\) |
| \((2,4)\) | 3 | \(5/8\) | 4 | \(0\) |
| \((2,4)\) | 4 | \(2/3\) | 3 | \(1\) |
| \((2,5)\) | 4 | \(3/5\) | 5 | \(0\) |
| \((2,5)\) | 5 | \(5/8\) | 4 | \(1\) |
| \((2,5)\) | 10 | \(2/3\) | 3 | \(2\) |

For reproducibility the table uses this rational enumeration.

~~~python
from fractions import Fraction

for ell, mu in [(4, 2), (4, 3), (2, 4), (2, 5)]:
    beta = mu * (ell - 1) - ell
    for degree in range(2, beta * mu + 2):
        s = 1 - Fraction(beta * degree, ell * mu * (degree - 1))
        if not Fraction(1, ell) < s < Fraction(1, ell) + Fraction(1, mu):
            continue
        h = (ell * s).denominator
        lowest = [j for j in range(mu) if (mu - j) % h == 0]
        if lowest:
            print(ell, mu, degree, s, h, lowest)
~~~

## Removing the remaining successors

For \((\ell,\mu)=(2,4)\), the row \(M=3\) gives a binomial
\(f=\eta^4-u\) with four distinct roots, too many for \(W\).
The row \(M=4\) gives \(f=\eta(\eta^3-u)\), with simple
roots and multiplier \((1-s)/w=M/\mu=1\) in (2), also
impossible. For \((2,5)\), the rows \(M=4\) and \(M=5\)
fail for the same reasons with one additional root.
In the row \(M=10\), \(f=\eta^2(\eta^3-u)\) with \(u\ne0\),
and the multiplier in (2) is \(2\).

Project \(W\) to its exponents congruent to one modulo three.
It still solves (2) with the same constant, because the
operator in (2) lowers that character by one under
\(\eta\mapsto\zeta\eta\), \(\zeta^3=1\). More explicitly
\(W(\zeta\eta)\), as a function of \(\eta\), has right
side \(\zeta c_0\), so average with coefficients \(\zeta^{-1}\).
This polynomial has degree ten, vanishes on the roots of
\(f\), and has the form \(W=\eta(t-u)V(t)\), where
\(t=\eta^3\) and \(V(t)=A+Bt+Ct^2\). Substitution gives
\[
                 t(t-u)V'+(u-2t)V=c_0/3 .
\]
The constant term \(uA\) must be nonzero, while the coefficient
of \(t\) is \(-2A\), a contradiction. This eliminates
the entire \((6,10)\) corner.

For \((4,3)\), the row \(M=3\) gives
\(f=\eta(\eta^2-u)\) and multiplier one.
The row \(M=6\) has multiplier two, so (2) would imply
\((W/f^2)'=c_0/f^2\). At any simple root \(\alpha\)
of \(f\), the residue of \(1/f^2\) is
\(-f''(\alpha)/f'(\alpha)^3\), which would have to vanish.
A cubic with three simple roots cannot satisfy this because
its second derivative has degree one. A cubic with a
double root \(a\) and a simple root \(b_0\ne a\) cannot
either, since in its monic form \(f''(b_0)=4(b_0-a)\ne0\).
Thus no successor with multiple distinct roots remains
in the \((4,3)\) case. A single-root successor at \(s=1/2\)
only returns the same \(C\) and denominator, and
Proposition 4.3 requires a strictly larger next slope.
It cannot avoid the exclusions just proved.

For \((4,2)\), the row \(M=2\) cannot have multiple distinct
roots: they would both be simple with multiplier one.
Thus, after possibly one single-root shift at \(s=1/2\),
the only remaining edge has \(s=5/8\), \(M=3\),
and \(f=\eta^2-u\) with \(u\ne0\).
Here \(C=(3/2,2)\) and \(w=1/4\).
At this edge Proposition 4.3 applies again by the
hypotheses checked above, with a simple root of \(f\)
of multiplicity two in the face of \(P\).
The denominator becomes \(\operatorname{lcm}(4,8)=8\)
and the new normalized end is
\[
                       C'=(1/4+5/8,\,1)=(7/8,1).
\]
The successor has slope \(s'>5/8\) with \(s'<7/8\)
by positivity, and ratio \(2:b\).
If its leading bracket is scalar, then
\((b+2)(7/8-s')=1-s'\), or
\[
                        8s'=7-\frac1{b+1}.
\]
The denominator \(b+1\ge4\) is larger than the maximal
\(Y\)-exponent \(2\) of the first face, so this face
would be a monomial. If the leading bracket is zero,
its common factor has \(Y\)-degree one and must have
two terms because it is an edge. This requires \(8s'\)
integral. The sole possible value is \(s'=3/4\).
Its nonzero root shift keeps \(C'\) and denominator
eight. All the hypotheses of Proposition 4.3 remain
valid, including the negative \((1,-1)\) endpoint weight
from \(7/8-1<0\). A further successor must have strictly
larger slope, and neither a scalar nor a commuting
successor is then possible. This final contradiction
eliminates the \((4,12)\) corner and proves the theorem.
