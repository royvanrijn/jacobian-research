# Stationary contacts unfolded by shifting a retained midpoint

Consider a simple interior Gaussian contact of an ordinary finite
or countably infinite retained core, with the conventions of the
[compactified interior reference](MULTIVARIATE_COMPACTIFIED_INTERIOR_CONTACT_FINITE_REDUCTION.md).
Assume the angular difference has even order,
\(\theta(T_*+x)=c x^p+O(x^{p+1})\), with \(c\ne0\), and both
individual angular velocities vanish at \(T_*\).

Let \(\mathfrak L_T(P)\) denote the first variation of the
angular difference when the coefficient at \(k\) of order \(j\)
is multiplied by \(\exp(v P_j(k/j))\) and \(v\) varies at zero.
Use the zero-endpoint shapes \(\phi,\psi,\Xi\) of the
[third-order chord formula](MULTIVARIATE_THIRD_ORDER_BOUNDARY_AND_INTERIOR_REDUCTION.md#exact-row-expansions-and-the-projected-coefficients),
and assume \(\mathfrak L_{T_*}(j\psi)=0\). Define the two
independent responses
\[
                  C=\mathfrak L_{T_*}(j\Xi),\qquad
                  B=\mathfrak L_{T_*}(\psi).
\]

**Theorem.** A sequence of simple sign-changing interior
contacts of finite ordinary cores tends to the given contact
in each of the following cases:

- \(cC<|cB|\).
- \(R\) is infinite and \(cC>|cB|\).
- \(R\) has infinitely many rows with signed gap different
  from \(\epsilon\), where \(cC=|cB|>0\) and \(C=\epsilon B\)
  with \(\epsilon\in\{1,-1\}\).
- There are infinitely many retained rows of one fixed signed
  gap \(d\), with \(dC=B\). Put \(v=s-1/2\) and
  \(\Upsilon(s)=v/320-v^5/20\).
  The fourth coefficient
  \[
       J_4^{(d)}=
         \mathfrak L_{T_*}\!\left(
               j d^5\Upsilon-\tfrac32 d^4\Xi+
                                      \tfrac{d^3}{2j}\psi\right)
  \]
  is nonzero.

These include all nonzero third-order response cases \(C,B\)
for an infinite retained core, except for the equality case
with eventually only the indicated aligned unit gaps. The
fourth coefficient test applies also in that equality case
and when \(C=B=0\) and a fixed gap repeats infinitely.
The consecutive orders and total integer degree difference
are preserved, along with every strict initial curvature cap.

## The two third-order responses

For any weak ordinary row with midpoint \(M\) and signed gap
\(d\), use the [uniform compensated log expansion](MULTIVARIATE_REVERSE_CLOCK_ROW_REPLACEMENT.md#uniform-compensated-row-expansion).
An individual multiplier proportional to \(\phi\) is a
change of Gaussian clock in that order, so its linear angular
response vanishes at \(T_*\). The signed leading cubic
response vanishes by hypothesis. Thus the remaining phase
change at \(T_*\) on adding or removing the row and
compensating its variance is
\[
          \pm\frac{d^4}{M^3}\left(C-\frac{B}{d}\right)
                  +o(d^4/M^3).                              \tag{1}
\]
The estimate holds also when \(|d|\to\infty\). The linear log
remainder is \(O_j(|d|^5/M^4)\), and the nonlinear root
remainder is \(O_j((|d|^3/M^2)^2)\). Their ratios to
\(d^4/M^3\) tend to zero with \(|d|/M\) and \(d^2/M\),
respectively.

In the first case choose a signed unit gap \(d\) for which
\(c(C-B/d)<0\) and add an ordinary row with that gap at a large
midpoint. In the second case delete a weak retained row
since \(c(C-B/d)\ge cC-|cB|>0\) for every nonzero integer \(d\).
In the third case delete weak retained rows with \(d\ne
\epsilon\), for which \(c(C-B/d)\ge cC/2>0\).
Infinitely many eligible rows always supply a sequence with
\(d^2/M\to0\) by summability. In each case (1) gives phase
strictly opposite to \(c\) at the original contact strength.

## Shifting one fixed-gap row

For fixed \(d\ne0\) write \(u=1/M\). The exact normalized
compensated log multiplier of a row is analytic in \(u\)
near zero. This follows by cancelling the common power of
\(M\) in the usual finite products for the nonconstant
coefficient ratios. Its first nonzero terms are of order
\(u^2,u^3,u^4\). The order-four coefficient is
\[
                     j d^5\Upsilon-\tfrac32d^4\Xi+
                                         \tfrac{d^3}{2j}\psi .
\]
Indeed in \(d^2j\psi_1(1+j(M+d v))\) the \(u^4\) term is
\(-d^5v^3-(3d^4/(2j))v^2-(d^3/(2j^2))v\).
The zero-endpoint primitives in the multiplier are
\(\Upsilon''=-v^3,\ \Xi''=v^2,\ \psi''=-v\), and there
is the overall factor \(j\) converting the chord to the log
coefficient.

Write \(J_h^{(d)}(T)\) for the linear phase response to the
coefficient of \(u^h\) in this multiplier. At \(T_*\),
\(J_2^{(d)}=0\) by the velocity and cubic assumptions,
and \(J_3^{(d)}=d^4C-d^3B\). Under the fourth hypothesis
both are zero. Move a weak retained midpoint from \(M\)
to \(M+\sigma\) for \(\sigma=1\) or \(-1\), using the same
gap and compensating the leading variance difference.
Both endpoints still have integer exponents. Since
\[
         (M+\sigma)^{-1}=u-\sigma u^2+O(u^3),
\]
the difference between the two compensated row multipliers
has norm \(O(u^3)\). The resulting phase difference at
\(T_*\) is
\[
                   -4\sigma J_4^{(d)}u^5+O(u^6).              \tag{2}
\]
In (2) the nonlinear angular remainder is \(O(u^6)\)
by analyticity of the simple-root map. The sign choice of
\(\sigma\) therefore gives phase strictly opposite to \(c\).

In all cases use finite truncations with coefficient errors
smaller than the displayed perturbations, including any
selected row before deleting or shifting it. The two
endpoint signs at nearby strengths around the original
contact persist under sufficiently small perturbation,
while the new phase at \(T_*\) has the opposite sign.
Continuity and real analyticity give a sign-changing zero
on a shrinking strength interval. The original root
simplicity persists. Recover the physical bias by matching
the positive root modulus ratio, keeping \(\Delta_0\) fixed
as in the [finite truncation proof](MULTIVARIATE_COMPACTIFIED_INTERIOR_CONTACT_FINITE_REDUCTION.md#recovering-the-physical-scale).
The row and center changes tend to zero, so the new contact
strength and bias converge and their strict interior
inequality persists.

## A transverse test by moving the same row

Suppose there are infinitely many rows of a fixed signed
gap \(d\), the angular order is even \(p\ge4\), and
\(\partial_T J_2^{(d)}(T_*)=b\ne0\). When \(p=4\), require
also \(J_3^{(d)}(T_*)=0\). The unit midpoint move then
gives a sequence of **transverse** interior contacts.
Indeed, uniformly for \(T\) near \(T_*\), the exact
phase after the move has expansion
\[
 \theta(T)-2\sigma u^3 J_2^{(d)}(T)
        +u^4\bigl(3\sigma^2 J_2^{(d)}(T)-3\sigma J_3^{(d)}(T)\bigr)
        +O(u^5)
\]
with one \(T\) derivative. At
\(T=T_*+u^{3/(p-1)}y\), divide by \(u^{3p/(p-1)}\).
The limit is \(c y^p-2\sigma b y\), including its first
\(y\) derivative. The \(u^4\) terms are negligible because
the exponent is below four for \(p>4\), and for \(p=4\)
both coefficients of \(u^4\) vanish at \(T_*\).
For a suitable sign \(\sigma\) the limiting equation has
a nonzero simple real root. Use truncation errors
\(o(u^{3p/(p-1)})\) and recover the physical scale as above.
