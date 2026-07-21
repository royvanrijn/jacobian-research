# Boundary-intersection obstruction for C24

This note settles the remaining comparison between `C24_(m,1)`, `m>1`, and
a generic weighted seed of inverse degree `m+2`.  The obstruction is the
reduced intersection of the two intrinsically distinguished target boundary
components.  The general theorem making this construction a canonical stable
left--right invariant is recorded in
[CANONICAL_BOUNDARY_INTERSECTION_INVARIANT.md](CANONICAL_BOUNDARY_INTERSECTION_INVARIANT.md).

Throughout, `k` is a characteristic-zero field.  For the weighted family,
write

\[
 H(W)=W^2J(W),\qquad J(0)=h_2\ne0,
\]

where `J` is squarefree, has no zero at the origin, and has the distinguished
simple zero at `W=1`.  These are the generic weighted seeds considered in the
resolvent--ramification signature.

## The two intersections

For C24 with `r=1`, the inverse polynomial is

\[
 \Psi(T)=C\int_0^T\{1-t(Q-Pt)^m\}\,dt-R,
 \qquad C\in k^*.
\]

Its discriminant divisor is the image of

\[
 1-T(Q-PT)^m=0.
\]

This critical divisor has coordinate ring `k[Q,Y,Y^(-1)]`, with
`T=Y^(-m)` and `P=(Q-Y)Y^m`; it maps birationally to the discriminant.
Its affine chart sees `Q=Y` when `P=0`; the target closure also retains the
limit cluster at `Y=0`.

On the second nonproperness component `P=0`, the finite critical chart gives
`T=Q^{-m}`.  Substitution into the critical value gives

\[
 R=C\left(T-\frac{Q^mT^2}{2}\right)=\frac{C}{2Q^m}.
\]

This is one component of the target intersection.  The closure of the
critical divisor contributes a second component at `Q=0`.  Uniformly, direct
expansion of the resultant defining the discriminant gives

\[
 \left.P^{-m(m-1)}\operatorname{Disc}_T(\Psi)\right|_{P=0}
 =\lambda C^{2m+1}Q^{m(m+1)}(2RQ^m-C),                \tag{1}
\]

where `lambda in k^*`.  Here is the resultant calculation of the lowest
part.  On the critical divisor put

\[
 F(Y)=P-(Q-Y)Y^m
\]

and multiply the critical-value equation by `Y^(2m)`.  It becomes

\[
 G(Y)=C\{Y^m-S(Y,Q)\}-RY^{2m},
 \quad
 S(Y,Q)=\int_0^1u\{Q+(Y-Q)u\}^mdu.
\]

At `P=0`, one has `F=Y^m(Y-Q)`, while

\[
 G(0)=-\frac{CQ^m}{(m+1)(m+2)},\qquad
 G(Q)=Q^m\left(\frac C2-RQ^m\right).
\]

Therefore `Res_Y(F,G)` is a nonzero scalar times
`C^m Q^(m(m+1))(C-2RQ^m)`.  Comparing the standard
discriminant--derivative resultant before and after `Y=Q-PT` gives the
additional nonzero power of `C` in (1); the leading coefficient of `Psi`, of
`P`-order `m`, gives the exact saturation order `m(m-1)`.  This proves (1),
not merely its reduced support.  Thus the reduced target intersection is

\[
 \boxed{P=0,\quad Q(2RQ^m-C)=0.}                       \tag{2}
\]

The two factors in (2) are comaximal because `C` is nonzero.  Its coordinate
ring is therefore

\[
 k[Q,R]/(Q(2RQ^m-C))
 \simeq k[R]\times k[Q,Q^{-1}],                        \tag{3}
\]

so the intersection is `A^1 disjoint-union G_m`.

For the weighted seed, the inverse polynomial is

\[
 E(W)=H(W)-BCW+cAC^2,\qquad c\in k^*.
\]

Put `W=Cu` near the double zero of `H`.  Dividing by `C^2` gives

\[
 h_2u^2-Bu+cA+O(C).
\]

All roots of `J` are simple and nonzero, so every discriminant factor
involving one of them is a unit at the generic point of `C=0`.  The two small
roots alone contribute the exact saturated trace

\[
 \boxed{C=0,\quad B^2-4h_2cA=0}                        \tag{4}
\]

up to a nonzero scalar.  Thus this reduced intersection has coordinate ring

\[
 k[A,B]/(B^2-4h_2cA)\simeq k[B].                       \tag{5}
\]

This root-scaling argument also proves that the pulled-back discriminant has
exact `C`-adic order two.  It is the uniform version of the one-extra-root
calculation in `verify_deformed_seed_boundary.py`.

## Intrinsic character and inequivalence

In both families the target discriminant component is intrinsically singled
out by its boundary ramification index two and sheet loss two.  The second
target component is intrinsically singled out by its unramified boundary
primes and total sheet loss `m-1>0`.  A polynomial left--right equivalence
must therefore carry these two components, respectively, to their
counterparts.  Its target automorphism would restrict to an isomorphism of
their reduced intersections.

No such isomorphism exists.  The ring (3) has a nontrivial idempotent, whereas
the domain (5) does not.  Equivalently, the C24 intersection is disconnected
and the weighted intersection is connected.

### Theorem

For every `m>1`, no `C24_(m,1)` map, including any allowed tail
`h_q(A)+A^2g(A)`, is polynomially left--right equivalent to a generic
weighted-seed map of inverse degree `m+2`.

The same conclusion holds after adjoining equally many identity coordinates,
by the canonical boundary-intersection theorem.  Concretely, the
distinguished intersections become
`(A^1 disjoint-union G_m) x A^a` and `A^(a+1)`; polynomial extension does not
remove the nontrivial idempotent.

The tail does not enter the inverse polynomial, its critical divisor, or
either target boundary component, so the argument is uniform in `g`.

Run the bounded exact regression with

```bash
.venv/bin/python scripts/verify_boundary_intersection_obstruction.py
```

It checks inverse degrees four through six.  The all-degree result is the
argument above, not an extrapolation from those cases.
