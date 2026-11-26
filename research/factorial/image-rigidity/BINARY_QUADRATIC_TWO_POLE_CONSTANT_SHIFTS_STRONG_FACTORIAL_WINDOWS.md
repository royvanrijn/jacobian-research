# Constant shifts of binary quadratics from residue poles and coefficient grids

Let K be a finite extension of Q2, with valuation v(2)=1, valuation ring O
and residue field F. More generally the proofs apply to a real-valued
nonarchimedean extension of v2. Let

    q=(alpha/2)X^2+beta XY+(gamma/2)Y^2+dX+eY,
    alpha,beta,gamma,d,e in O,
    R_j=Lambda(q^j)/j!,  R_0=1.

The [noncritical companion](QUADRATIC_CONSTANT_NONCRITICAL_DYADIC_WINDOW_TRANSFER.md)
proves these coefficients integral, and gives their universal formal residue

    Rbar(T)=N(T)/D(T),
    N=1+(dbar+ebar+betabar)T+dbar ebar T^2,
    D=(1+alphabar T+dbar^2 T^2)(1+gammabar T+ebar^2 T^2).       (1)

All poles below mean roots of the REDUCED denominator over an algebraic
closure of F. In particular they are finite geometric poles, not a claim
about analytic singularities. No Hessian-rank or fixed-correlation equation
is assumed. All powers and generating series are formal.

**Theorem 1 (uniform two-pole constant-shift six-window).** If Rbar has at
least two distinct geometric poles, then for EVERY k in K, including k=0,
the phase f=q+k has no six consecutive zero factorial moments at any
positive start. The constant k need not be integral and no fixed constant
valuation or field ramification index is required.

For k=0 there is in fact a four-window. For nonzero k this theorem does
not assign a phase with fewer than six monomials its smaller actual-support
window merely from the six-window.

**Corollary 2 (dense five/six-support SFC sector).** Suppose all five
coefficients alpha,beta,gamma,d,e are nonzero and Rbar has at least two
distinct geometric poles. Then q+k satisfies its ACTUAL-support strong
factorial window for every k in K: the actual support is five when k=0
and six when k!=0. This is a local arithmetic sufficient condition for
that full dense family, not a proof of arbitrary binary-quadratic SFC or
SFC(2). There is no assertion that every complex phase has an embedding
and normalization with two residue poles.

## 1. Proof over every constant valuation

Two poles imply that Rbar is nonpolynomial. Section5 of the noncritical
companion proves a positive unit-coefficient window of length
w=deg(reduced denominator)<=4, including the equal numerator/denominator
degree case. For k=0, M_m=Lambda(q^m)=m!R_m is nonzero whenever R_m is a
unit. Hence the four-window, and in the dense case the five-support
window, follow immediately.

Assume k!=0 and put delta=v(k).

* If delta<0, every positive moment is nonzero.
* If delta=0, every positive even-order moment is nonzero, giving a
  two-window. This does not assert transfer of a one-window.
* If 0<delta<1 is not a critical number1-2^(-a), the unique minimizing
  factorial index and circular-prefix transfer in the companion give
  the same w-window, in particular a four-window.
* If delta>=1, the companion's unit-window transfer again gives the
  same w-window.
* At the first critical value delta=1/2, its separate elementary
  multiple-of-four argument gives a four-window for every integral
  divided-factorial sequence, independently of the pole hypothesis.
* At any remaining critical value delta=1-1/L, L=2^a,a>=2, let
  mu be the nonzero residue of k^L/2^(L-1). The polynomial
  1+mu T^L has just one geometric root. It cannot cancel two distinct
  poles of Rbar, so B=(1+mu T^L)Rbar is nonpolynomial. The
  [critical companion](BINARY_QUADRATIC_CRITICAL_DYADIC_SIX_MOMENT_WINDOWS.md)
  proves the six-window, using the exact native numerator relation.

These cases exhaust every real-valued delta. When 0<delta<1, this uses
the elementary partition into the open bands
1-2/2^a<delta<1-1/2^a and their endpoints; no higher-rank ordered-value-group
exhaustion is asserted. This proves Theorem1. In the dense case the k=0
four-window fits within the actual five-window and each k!=0 six-window
is its actual support-size window, proving Corollary2.

## 2. A sharper per-constant statement and the exact uncovered frontier

The two-pole condition is a convenient UNIFORM sufficient condition over
all constants. For a specified nonzero constant in a noncritical open
band or with delta>=1, mere nonpolynomiality of Rbar suffices, and the
companion gives a window no larger than the actual support of q+k.
At a specified critical value with a>=2 it suffices that

                (1+mu T^L)Rbar(T) is nonpolynomial.          (2)

This includes one-pole residues whose pole is not completely canceled.
No converse from failure of these criteria to an actual moment zero block
is asserted.

**Corollary 3 (necessary residual exception for a six-zero block).** Under
the integral natural-coefficient setup, if q+k has a six-zero block then
either Rbar is polynomial, or k!=0, delta=1-1/L with L=2^a,a>=2, and
the reduced denominator has a single geometric pole completely canceled
by 1+mu T^L. In the nonpolynomial exceptional case it can be written

            D0(T)=(1+rho T)^h,  1<=h<=4<=L,
            mu=rho^L                                      (3)

over the algebraic closure of F. For k!=0 a six-zero block also requires
delta>=3/4. This is a NECESSARY exception list, not a list of counterexamples.
It does not remove the smaller-support critical-window question.

For the proof, if Rbar is nonpolynomial and k=0, the four-window already
excludes the block. The same is true at every noncritical positive
valuation and at delta>=1. Negative, zero and the first half-critical
valuation are excluded without any residue-pole assumption. In the open
band1/2<delta<3/4, the noncritical theorem has L=4 and gives a four-window
from R_0=1 even without a residue unit-window hypothesis. Hence that band
is also excluded, so a six-zero block with k!=0 requires delta>=3/4.
At a higher critical radius the critical theorem excludes a block unless B is
polynomial. In a reduced fraction N0/D0, N0 and D0 are coprime, so B
polynomial forces D0 to divide1+mu T^L. Since this last polynomial has
one geometric root, D0 has one pole and multiplicity h<=degD<=4.
As L>=4, its entire denominator is canceled exactly in the situation(3).
This proves the stated necessary condition.

## 3. Coefficient value grids: a second all-constant sector

The two-pole hypothesis in Theorem1 is independent of a separate coefficient-grid sufficient condition. In this section let E be ANY characteristic-zero nonarchimedean valued field, v(2)=1, whose coefficient value group satisfies

                  v(E*) contained in (1/e0)Z               (4)

for a positive integer e0. Take alpha,beta,gamma,d,e in O_E and let the constant k lie in ANY valued extension of E. The natural nonconstant coefficients, and therefore every integral R_j, belong to E; the constant need not belong to E or its value grid. If a moment vanishes, k is algebraic over E and a least-term tie in its minimal polynomial puts v(k) in the rational divisible hull of (4), even for a higher-rank extension valuation.

**Theorem 4 (small-grid all-constant windows).** Suppose Rbar is nonpolynomial. If e0 is odd, q+k has no five consecutive zero factorial moments. If e0 is twice an odd integer, it has no six consecutive zero factorial moments. These statements hold for every constant k in every valued extension of E. In particular, 4 does not divide e0 is sufficient for a six-window, without any pole-count or Hessian-rank hypothesis.

The [multiplicity companion](QUADRATIC_CONSTANT_APPELL_MULTIPLICITIES_FROM_DISCRETE_VALUE_GRIDS.md) proves the more general transfer from a positive unit window w to w+2^v2(e0). Its argument is finite: the least-cost exponent width at a critical radius theta=1-1/L is exactly L whenever its selected anchor is a unit and a root exists. Widths add under polynomial multiplication. The coefficient grid forces the positive width of the minimal factor to be at least L/gcd(L,e0), so the root multiplicity is at most gcd(L,e0). Appell differentiation excludes a w+gcd(L,e0) zero block. The native binary residue has w<=4 by the coefficient companion, proving Theorem4. No Newton root-count or completeness theorem is used.

When all five nonconstant coefficients are nonzero, the k=0 four-window and the k!=0 six-window prove the ACTUAL-support SFC for this dense family under either sufficient hypothesis: two geometric residue poles as in Corollary2, or nonpolynomial residue and 4 not dividing e0 as in Theorem4. The odd-grid five-window is an additional stronger uniform bound; sparse phases are not assigned their smaller actual-support bounds merely by deleting monomials.

**Corollary 5 (coefficient-grid obstruction to a counterexample).** In the setup(4), if q+k has six consecutive zero factorial moments and Rbar is nonpolynomial, then 4 divides e0. Moreover k!=0, delta=v(k)>=3/4 is a higher critical value 1-2^(-a), a>=2, and the one-pole complete-cancellation conditions(3) hold. This is a necessary condition, not an existence assertion.

The first conclusion is the contrapositive of Theorem4. A possible zero makes delta rational, as explained above, so the proof of Corollary3 applies with exactly the same native moment identities and gives every other condition. For a given coefficient field the assertion applies to each admissible grid bound e0. It does not assert that every phase has a normalization over an odd or twice-odd grid, and it does not identify an arbitrary grid bound with the least ramification index. Polynomial residue series remain outside Theorem4 and outside this corollary's nonpolynomial hypothesis.

## 4. Normalization and scope

For a nonconstant binary quadratic Q over a chosen 2-adic coefficient
field, divide the WHOLE phase Q+k by an element realizing the minimum
of the valuations of its five natural nonconstant coefficients
2A,B,2C,d,e. The resulting alpha,beta,gamma,d,e are integral and at least
one is a unit. The constant is divided by that same element. Whole-phase
scaling preserves every moment zero and every nonzero monomial support.
The pole test and relative constant valuation are computed after this
normalization. A constant-only phase has no zero positive moments and
requires no such normalization.

The proofs use only the complete coefficient identities in the three
companions, factorial valuations, multiplicative least-cost exponent widths,
finite leading-residue calculations, and recurrences for formal rational series. They do not import a
matching-family upper theorem, a Gaussian holonomic-rank window, a
Newton-factor precision conjecture, or a finite experimental extrapolation.
The polynomial-residue and complete one-pole-cancellation exceptions
remain explicit. All conclusions concern the specified valued sectors;
the unrestricted all-start binary-quadratic, SFC(2) and SFC(3) goals are
not closed by them.
