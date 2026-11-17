# Field trace preserves divergence, not automatically the branch curve

14 November 2026. PF2TRACE1 proves the valid trace consequence of the
Wright grading and records an exact finite quartic control rejecting a
stronger, unjustified inference. The trace vector field is polynomial
with divergence four. Tangency to the branch curve is an additional gate;
it does not follow just from source tangency to the ramification prime,
semisimplicity, or homogeneity of the pulled-back volume.

## The valid finite-cover trace statement

Let pi:N -> A2_(p,q) be finite dominant from a normal affine complex
surface, with function fields K/F of degree d. Let delta be a regular
derivation of Gamma(N), and write sigma=dp wedge dq on the generic point.
Assume Lie_delta(sigma)=kappa*sigma for a constant kappa. Then

\[
 \Theta=\operatorname{Tr}_{K/F}(\delta p)\partial_p+
        \operatorname{Tr}_{K/F}(\delta q)\partial_q
\]

is a polynomial derivation of C[p,q] and has divergence d*kappa.
The coefficients are polynomial because delta(p),delta(q) are integral
over the target and their traces belong to its integrally closed fraction
field. No unproved degree bound on these polynomial coefficients is used.

Over a characteristic-zero separable field extension, partial derivatives
on F extend uniquely to K and commute with trace. For example they extend
uniquely to a Galois closure and commute with every F-automorphism, so
one may differentiate the sum over embeddings. Writing
 delta=(delta p)*partial_p+(delta q)*partial_q on K gives

\[
 \operatorname{div}(\Theta)
 =\operatorname{Tr}(\partial_p\delta p+\partial_q\delta q)
 =\operatorname{Tr}(\kappa)=d\kappa.
\]

For the actual Wright ring, Delta=-X*partial_X+2Y*partial_Y is regular,
preserves R, and has divergence one with respect to the prescribed
source volume. Thus any finite quartic realization supplies such a
polynomial Theta with divergence four. This alone is not a branch-
preserving Euler field.

## Where an attempted tangency argument loses the affine sheets

Let h=0 define the branch curve, whose generic fibre has the one ramified
prime R with e=2,f=1 and an unramified algebra of total rank two. The
reduction of the finite trace modulo h is the ramified contribution with
multiplicity two plus the unramified traces. Since delta is tangent to R,
its ramified contribution to the normal derivative of h is zero. Therefore

Theta(h) mod h = the trace of delta(h) on the unramified rank-two part.

There is no reason for this remaining trace to vanish. Source preservation
of R does not say that delta preserves the other primes over the same
branch divisor. Proving the unramified normal trace vanishes would be
an additional constraint; it is not supplied by constant divergence.

## A fully explicit finite quartic countercontrol

This control has N=A2_(x,w), not the Wright normalization. It is not a
Keller counterexample. Put

\[
 p=x+w+x^3,\qquad q=-x^2/2-xw-3x^4/4.
\]

Its Jacobian is w. The ramification divisor is the single smooth line
R=(w=0), generically of index two. The map is finite free of degree four:

\[
 \mathbb C[x,w]\simeq
 \mathbb C[p,q,x]/(x^4+2x^2-4px-4q),\qquad
 w=p-x-x^3.
\]

On R the branch parametrization is

p(t)=t+t^3, q(t)=-t^2/2-3t^4/4.

Its rational inverse is

 t=p(1-12q)/(1+9p^2+4q),

so the generic residue degree is one. Thus the control retains finite
degree four, one reduced ramification prime, its affine-line normalization
and the two generic unramified sheets.

Use the regular semisimple hyperbolic grading field

 delta=(-x*partial_x+2w*partial_w)/3.

It preserves R and has Lie_delta(w dx wedge dw)=w dx wedge dw.
The exact power traces of the monic quartic give

Tr(x)=0, Tr(x^2)=-4, Tr(x^3)=12p.

Since delta(p)=(2p-3x-5x^3)/3 and
 delta(q)=(16q-6x^2+15px)/3, its target trace is

\[
 \boxed{\Theta=-\frac{52}3p\partial_p+
            \left(\frac{64}3q+8\right)\partial_q.}\tag{1}
\]

The divergence is four. At t=0 the branch point is(0,0) and its tangent
is horizontal, whereas Theta(0,0)=(0,8). Thus Theta is not tangent to
the branch curve. The exact discriminant check gives the same conclusion.

At this target point the remaining roots are x^2=-2, with w=x. Each is
unramified and contributes four to the trace of delta(q); the ramified
point contributes zero. These are exactly the affine-sheet terms that
an invalid tangency proof would discard.

The countercontrol has a principal canonical/ramification class and its
complement of R is not an affine plane. It therefore does not settle
trace tangency for the actual Wright candidate. It demonstrates that the
listed trace, grading and local ramification hypotheses alone are
insufficient; no actual Wright tangency theorem is claimed.

## Correct next-use boundary

A branch Euler or cusp-filtration obstruction cannot be applied from this
trace construction until the missing tangency condition is proved.
Separately, the [grading obstruction](WRIGHT_GRADING_OBSTRUCTION.md)
shows that the original source grading cannot itself preserve the target
Keller algebra. If a finite quartic candidate exists, at least one of

Delta(P)-Tr(Delta(P))/4, Delta(Q)-Tr(Delta(Q))/4

is a nonzero trace-zero element. The polynomial trace is an average; it
must not be substituted for a descended copy of the source grading.

The [producer](compute_wright_grading.py),
[certificates](wright_grading_certificates.json) and
[checker](verify_wright_grading.py) retain the monic quartic multiplication
matrix, power traces, exact discriminant, branch normalization and
transversality. The universal trace proof is written and is not formally
verified or independently reimplemented as a whole. Actual finite quartic
maps on the Wright model, higher poles and general ruled compatibility
remain open.
