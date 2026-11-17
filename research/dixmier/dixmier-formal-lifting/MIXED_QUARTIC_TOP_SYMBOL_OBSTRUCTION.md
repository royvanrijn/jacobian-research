# A non-square mixed quartic top symbol excludes every finite mate

Let k have characteristic zero, eta in k^times, and

\[
 P=\sum_{j=0}^4p_j(V)A^j\in k(V)[A;\eta V^{-1}\partial_V].
\]

Assume rho,g in k^times and sigma in k, with

\[
 \begin{split}
 p_4&=\rho V^{-6}+\text{terms of degree less than }-6,\\
 p_3&=\sigma V^{-3}+\text{terms of degree less than }-3,\\
 p_2&=g+\text{terms of negative degree},\\
 \deg_\infty p_1&\le2,\qquad \deg_\infty p_0\le5,\\
                         \sigma^2-4\rho g&\ne0.
 \end{split}                                             \tag{1}
\]

Then the rational-coefficient finite centralizer is k[P], and no
finite rational-coefficient operator has nonzero scalar commutator
with P. Mate nonexistence extends to every finite algebraic coefficient
cover. The centralizer equality itself is asserted over k(V).
There is no partner-order or rational position-degree bound.

This extends [QOI1's strict quartic cone](ORDERED_QUARTIC_INFINITY_OBSTRUCTION.md)
to a top cubic term. The discriminant assumption is essential to
this proof. Perfect-square top symbols require a different argument.

## Rational leading coefficients and finite even-flow exhaustion

If [P,T] has order at most two and T has positive order M, its highest
commutator coefficient gives

\[
            4p_4\delta(t_M)-Mt_M\delta(p_4)=0,
                     \qquad\delta=V^{-1}\partial_V.      \tag{2}
\]

At infinity, p_4 has degree -6. Hence t_M has degree -3M/2, so
M must be even for a rational leading coefficient.
If p_4 is not a square modulo a scalar, every possible M is divisible
by four: for M=4l+2, equation (2) would make p_4 a scalar multiple
of (t_M/p_4^l)^2. At an order divisible by four, subtract a constant
multiple of a power of P. Iteration ends in a constant. Thus in this
case the centralizer is k[P] and a nonzero scalar mate is impossible
without needing fractional positive powers.

Otherwise, after a harmless algebraic extension of constants, choose
an actual rational square root of p_4. The formal pseudodifferential
square root R=P^{1/2}, with that leading coefficient and order two,
has rational position coefficients, determined recursively by division
by 2sqrt(p_4). Set B_{2n}=(R^n)_+, n>=1. Every commutator [P,B_{2n}]
has order at most two, by its negative-part expression. Repeated
leading subtraction using (2) proves the exact finite normal form

\[
                    T=\sum_n c_nB_{2n}+c_0,\qquad c_n\in k.
                                                                  \tag{3}
\]

At order zero, a nonconstant coefficient would have commutator order
three, so the final coefficient is constant. If n is even, B_{2n}=P^{n/2}
and its commutator is zero. Only odd n need exclusion.
The final polynomial-centralizer conclusion descends from the constant
extension: compare highest A coefficients in a polynomial in P, subtract
that term and repeat. Each coefficient is rational in V and constant,
so belongs to the original k.

## The top Poisson bracket cannot vanish

Give V weight one and A weight three. Each nontrivial Ore contraction
loses five units of weight. With t=xi/V^3, the top symbol of P is

\[
          p_*=V^6F(t),\qquad
                  F(t)=\rho t^4+\sigma t^3+g t^2.        \tag{4}
\]

Every omitted term in (1) has strictly smaller weight. Root recursion
and finite multiplication show that the top symbol of B_{2n} is

\[
          b_*=V^{3n}G_n(t),\qquad
                   G_n=(F^{n/2})_+                       \tag{5}
\]

where the positive part is taken in the Laurent expansion at t=infinity.
This polynomial is nonzero and has degree 2n. Equations (4)--(5)
are coefficientwise weighted statements: lower input terms have lower
weight, and contractions lower weight by five. No unbounded spatial
series is evaluated.

The first contraction gives the top possible commutator symbol

\[
       3\eta V^{3n+1}\bigl(nF'(t)G_n(t)-2F(t)G_n'(t)\bigr).
                                                                  \tag{6}
\]

For example differentiation at fixed xi gives
partial_V(V^{3n}G_n(xi/V^3))=3V^{3n-1}(nG_n-tG_n'),
which verifies the coefficient and sign in (6).
If its polynomial factor vanished, then (G_n^2/F^n)'=0, hence
G_n^2=aF^n for a nonzero scalar a over an algebraic closure of k.
But F=t^2(rho*t^2+sigma*t+g), and the quadratic factor has two
simple nonzero roots by (1). For odd n their multiplicities in F^n
are odd, so it cannot be a polynomial square. Therefore (6) is nonzero
for every odd n.

Its weight 3n+1 is positive. In any nonempty finite odd part of (3),
the largest n has a nonzero top commutator weight strictly greater
than every lower odd flow; successive odd n differ by at least two,
so the weights differ by at least six. A scalar commutator, including
zero, cannot have this term. All odd coefficients in (3) vanish.
The remaining powers of P commute, proving both assertions over k(V).
A normalized coefficient trace gives mate nonexistence on every finite
algebraic coefficient cover.

## A finite parameter-pole descent corollary

Suppose an actual finite position admits a finite quantum rescaling
whose integral reduction is a P satisfying (1). A hypothetical finite
rational-position mate has finitely many poles in the scalar parameter's
Gauss valuation. Scale its nonzero scalar commutator to a unit.
At a highest parameter pole, its residue centralizes P and therefore
is a finite polynomial q(P). Subtract the same polar multiple of q
of the integral rescaled position. This preserves the commutator and
strictly lowers the pole order. After finitely many steps an integral
mate would reduce to a nonzero scalar mate of P, just excluded.
The argument also permits an algebraic-cover mate, traced down first.

Every subtraction is finite. No generic order bound or infinite
polynomial in the position is used. This is the finite pole-removal
argument of [QCL1](CRITICAL_LAURENT_FACE_OBSTRUCTION.md), with the
centralizer theorem above supplying its exact hypothesis. It is a
conditional same-critical-face statement, not a theorem about all
exact-R global images.

## A quadratic lower perturbation also has a universal residue obstruction

The following separate family has perfect-square top symbol, and so is
outside (1)'s discriminant hypothesis. It nevertheless has an independent
all-order obstruction. Let beta,g,eta be nonzero constants, mu arbitrary,
and tau!=0, and put

\[
 q=A+\beta V^{-3}A^2,\qquad
 P_\tau=gq^2+\mu V+\tau V^{-5}A^2.                       \tag{7}
\]

Then Cent(P_tau)=k[P_tau] over rational position coefficients, and no
finite nonzero scalar mate exists, including on finite algebraic covers.
This is not the previously closed first-order perturbation tau*V^{-2}A.

Adjoin a constant a_0 with a_0^2=g, set s=a_0*beta, and give V weight
one and A weight minus two. The lowest square part is
L_-^2, L_-=sV^{-3}A^2, of weight -14. The extra tau term has weight
-9, a gap of five. The linear A part in a_0*q also has gap five;
its complete contribution at tau=mu=0 is the actual differential
root a_0*q, so it contributes no negative part. The mu term has gap
fifteen. Thus the leading commutator of each nonpolynomial even flow
is linear in tau and calculated around L_- alone. Any additional
A part, tau insertion or mu insertion has strictly greater weight.

With theta=V partial_V and c(x)=x(x-2), the first tau-response of the
square root has the form

\[
 R_\tau=\frac\tau sV^{-2}F(\theta),\qquad
 c(x-2)F(x)+c(x)F(x-7)=c(x),\quad
                       F(x)=\tfrac12+x^{-1}+O(x^{-2}).   \tag{8}
\]

Indeed L_-R_tau+R_tau L_-=tau*V^{-5}A^2. Define the formal coefficient
functional mathcal L(B)=[x^{-1}]F(x)B(x), and
p_n(x)=product_{ell=1}^n c(x+7ell). Translation invariance of the
formal residue gives

\[
 \mathcal L\{c(x-2)B(x)+c(x+7)B(x+7)\}=0.               \tag{9}
\]

The derivative of P_tau^{k+1/2} is the sum L_-^i R_tau L_-^j,
i+j=2k, at that least weight. Translating each residue and successively
using (9) gives

\[
                     D_k=\mathcal L(p_{2k}),\qquad D_0=1. \tag{10}
\]

More explicitly, with b_i=product_{ell=0}^{i-1}c(x-2-7ell), relation
(9) gives mathcal L(b_i*p_j)=(-1)^i mathcal L(p_{i+j}).
The alternating 2k+1 translated terms therefore sum to (10).

For symbolic n the exact polynomial certificate is

\[
\begin{split}
 c(x+7(n+1))c(x+7(n+2))-r_n
   &=c(x-2)Q_n(x)+c(x+7(n+1))Q_n(x+7),\\
 Q_n(x)&=\frac{x^2}{2}+\frac{7(3n+4)}2x
                            +\frac{7(n+1)(35n+58)}4,\\
 r_n&=-\frac{7(n+1)(7n+9)(7n+11)(7n+16)}4.
\end{split}                                               \tag{11}
\]

Multiplying by p_n makes the right side (9)'s expression with B=p_n Q_n.
Consequently

\[
                      D_k=r_{2k-2}D_{k-1}\ne0
                                      \quad(k\ge1).      \tag{12}
\]

Every factor is nonzero in characteristic zero. This establishes all
indices by induction, rather than extrapolating a finite prefix.

The coefficient of A^{-1} in the least tau-response of P_tau^{k+1/2}
is D_k*tau*s^{2k-1}*eta^{4k+1}*V^{-4-14k}.
Only that negative coefficient can contribute to A^2 in minus its
commutator with P_tau. Direct leading multiplication therefore gives

\[
 [P_\tau,(P_\tau^{k+1/2})_+]_{A^2}
   =(56k+22)D_k\tau s^{2k+1}\eta^{4k+2}V^{-12-14k}
                         +\text{terms of greater degree}.        \tag{13}
\]

The numerical factor is -[4(-4-14k)-6]. For k=0 the positive square
root is B_2=a_0*q+tau*V^{-2}/(2s), and exact Ore multiplication gives

\[
 [P_\tau,B_2]_{A^2}
             =-a_0\eta\tau(V^5-22\beta\eta)V^{-12},     \tag{14}
\]

which independently checks the nonzero leading coefficient in (13).
Identity (11) was also checked as a polynomial in x and n, not by
substituting a finite index range.

The leading coefficient of P_tau is s^2V^{-6}, so every rational
finite scalar-commutator candidate has even order. Leading subtraction
is exhaustive and gives a polynomial in P_tau plus finitely many
positive half-integer flows. Their distinct nonzero lowest coefficients
(13) cannot cancel. Thus every scalar commutator, including zero,
forces the half-integer-flow coefficients to vanish. This proves the
centralizer assertion and mate nonexistence; normalized coefficient
trace supplies the finite algebraic-cover assertion for mates.

For a concrete global application, the exact admitted perturbation
S_lambda+nu*h^2*X^5*C^2 of the original mixed global S_lambda has
critical reduction (7), with g=-3r, beta=lambda*r^2, mu=2 and tau=nu*r^2.
The addition is a finite global centralizer element and is h-divisible.
The critical reduction follows directly from h=t^5,v=t^2V,C=rt(A-tV).
For h-independent lambda,nu!=0, finite parameter-pole descent using
the centralizer assertion excludes every finite generic mate of this
specified global perturbation. This does not classify arbitrary lower
momentum terms on perfect-square top faces.

## Scope and remaining structural requirement

The first theorem closes the entire stated strict cone when the top
quartic is not a square. Appendix (7)--(14) closes one additional
perfect-square-top quadratic perturbation by an exact universal
moment certificate. Together with the existing square, adjoint and
first-order-perturbation theorems, these results leave a prospective
new face needing all of: perfect-square top geometry, a resonant
bottom divisor, and an essential nonselfadjoint lower momentum part
outside the proved moment-certificate families. No explicit new face
is asserted here to meet every remaining obligation.

This note contains written universal proofs and exact symbolic algebraic
identities, not a finite mate ansatz or a bounded search. It claims no
bounded Long lift, DC1/DC2 resolution, formal verification, independent
whole-proof replay, external review or novelty. No unconditional global
classification is inferred beyond the specific rescaling hypotheses
and the displayed admitted family.
