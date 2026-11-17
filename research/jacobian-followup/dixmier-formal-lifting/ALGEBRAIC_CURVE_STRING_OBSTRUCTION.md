# The order two and three Weyl relation on every algebraic curve

14 November 2026. **QSO1.** Let F be a one-variable algebraic function
field over a characteristic-zero constant field, let delta be a nonzero
derivation, and let F[D;delta] be the ring of scalar differential operators,
with Df-fD=delta(f). There are no operators P,Q of exact orders two and
three satisfying
\[
                         [P,Q]=\kappa,\qquad
                    0\ne\kappa\in\ker(\delta).       \tag{1}
\]

The leading coefficients may be arbitrary nonzero algebraic functions.
The coefficient curve, its genus, finite covers, poles and coefficient
degrees are unrestricted. The proof reduces (1) to an algebraic
Painleve-I solution and excludes it using local Puiseux exponents and
Riemann--Hurwitz. It is independent of prescribed arithmetic support.

[Kanel-Belov--Fedorov, Lemma 4.8](https://arxiv.org/html/2610.33919v1)
proves the monic rational case for delta=v^(-1) d/dv, on its specific
quadratic cover. The argument below removes that coefficient-field and
monicity restriction. No novelty priority is claimed for this classical
differential-operator/Painleve connection.

## Normalize the leading coefficients on a finite cover

It suffices to extend the constants to an algebraic closure C and work
in a one-variable function field over C. In characteristic zero a
derivation extends uniquely to algebraic extensions. Since the derivation
is nonzero, its constant field in a one-variable function field over
algebraically closed C is exactly C.

Write the highest terms of P,Q as a D^2 and b D^3. The coefficient of
D^4 in their commutator is
\[
                           2a\,\delta b-3b\,\delta a.
\]
Equation (1) forces this to vanish, so b^2/a^3 is constant. Adjoin a
square root rho of a and set E=rho D. The new coefficient derivation is
epsilon=rho delta. In this Ore presentation P is monic of order two,
and the leading coefficient of Q is b/rho^3, a nonzero constant.
Divide Q by that constant. The commutator remains a nonzero constant.

The Ore automorphism E -> E-ell_1/2, fixing the coefficient field,
removes the linear term of P. We may therefore assume, on a finite
algebraic cover, that
\[
 L=E^2+u,\qquad Q=E^3+b_2E^2+b_1E+b_0,\qquad
                              [L,Q]=\kappa\ne0.     \tag{2}
\]
Primes below mean epsilon. A change of algebraic curve or nonmonic
leading coefficients cannot evade this normalization.

## The commutator creates an algebraic time coordinate

Comparing coefficients in (2) gives b_2'=0. Subtract b_2 L from Q.
The next two coefficients then give
\[
 Q=E^3+\left(\frac32u+c_1\right)E+\frac34u'+c_0,
                         \qquad c_0,c_1\in C.
\]
The scalar commutator is
\[
                  -\frac14u'''-\left(\frac32u+c_1\right)u'
                               =\kappa.             \tag{3}
\]
Set w=u+2c_1/3. Then
\[
                        (w''+3w^2)'=-4\kappa.
\]
Crucially, no pre-existing primitive of 1 is needed. The relation itself
defines one inside the algebraic coefficient field:
\[
                  \tau=-\frac{w''+3w^2}{4\kappa},
                         \qquad \epsilon\tau=1.     \tag{4}
\]
Thus tau is nonconstant, and the coefficient field is finite algebraic
over C(tau). With y=-w/2, equations (3)--(4) become
\[
            \frac{d^2y}{d\tau^2}=6y^2+2\kappa\tau.   \tag{5}
\]
The function y would be algebraic over C(tau).

The automatic time coordinate is what permits arbitrary derivations:
even if epsilon originally had no evident algebraic slice, (1) would
supply it.

## No algebraic solution of the resulting Painleve equation

We prove the slightly more general statement
\[
                    y''=6y^2+\alpha\tau+\beta,
                         \qquad\alpha\ne0           \tag{6}
\]
has no algebraic solution over C(tau). This supplies a self-contained
algebraic proof of the classical nonexistence fact.

At any finite tau_0 put z=tau-tau_0 and take a Puiseux expansion of an
algebraic candidate y.

If y has a pole, balancing the lowest exponents of y'' and 6y^2 forces
\[
                             y=z^{-2}+\cdots.
\]
Indeed the pole exponent must be -2 and its coefficient must be one.
If a first nonintegral exponent mu occurs, then mu>-2. Its lowest
contribution to (6) has coefficient
\[
              a\bigl(\mu(\mu-1)-12\bigr)z^{\mu-2}.
\]
All earlier exponents are integers, so this term cannot be supplied
by them. The indicial equation factors as
\[
                       (\mu-4)(\mu+3)=0.
\]
Both roots are integers, contrary to the choice of mu.

If y is finite at tau_0, a first nonintegral exponent mu>0 contributes
a mu(mu-1) z^(mu-2) on the left. All nonlinear contributions with that
fractional part have strictly larger exponent. Hence mu=0 or 1,
again impossible for a nonintegral exponent.

Every local Puiseux expansion at a finite point therefore has only
integer powers. Apply this to the smooth projective curve of the minimal
field C(tau,y). Its map to P1_tau is unramified above every finite
point: y generates the extension, and all its local branches lie in
C((tau-tau_0)).

Let d be the degree of that map, g its genus, and r the number of
points above infinity. All ramification lies above infinity, with
total ramification d-r. By
[Riemann--Hurwitz](https://stacks.math.columbia.edu/tag/0C1B),
\[
                         2g-2=-2d+(d-r)=-d-r.
\]
Since d>=1, r>=1 and g>=0, this forces d=r=1 and g=0. Thus y is
rational in tau.

Finally, a rational y with a pole of integer order m>=1 at infinity
makes y^2 have degree 2m, greater than both the degree m-2 of y''
and the degree one forcing term in (6). If y is bounded at infinity,
the nonzero linear forcing term dominates. Both cases are impossible.
This proves the assertion and contradicts (5).

## Scope and useful sharpness checks

The result survives arbitrary finite algebraic extensions of the
coefficient field, changes of the curve parameter, algebraic leading
coefficients, and arbitrary rational pole divisors. It is not a finite
coefficient-degree test or a result requiring a specific quadratic cover.

The orders two and three are substantive hypotheses. In C(z)[D],
[D,z]=1, the pair
\[
                    P=D^2+z,\qquad Q=-D
\]
has commutator one and orders (2,1). Replacing Q by Q+P^2 gives orders
(2,4) and keeps the same commutator. Thus no blanket claim excluding
all pairs with positive differential orders follows.

Matrices and infinite differential-operator series are outside the
scalar finite-order statement. The right side of (1) must be a
nonzero differential constant. An h-adic coefficient series in
K(v)[[h]] need not be algebraic over K((h))(v), so even a fixed
formal operator-order bound alone does not place it under this theorem.
This result is not an unrestricted formal obstruction of the kind
already excluded by QFL1.

## Application to the Dixmier construction gap

Any proposed finite quantization that reduces on an algebraic curve
to a scalar canonical pair of exact differential orders two and three
is impossible. Altering coefficient degrees, allowing more poles,
passing to a higher-genus coefficient curve, or adjoining finite
algebraic functions cannot repair that route.

For Long's classical map, the localized scalar fiber profile in the
recent nonlifting paper has orders two and three. QSO1 excludes every
finite quantum realization retaining such a scalar reduction, without
having to prescribe a Frobenius center map. This does not prove that
every bounded full lift retains the profile or admits that reduction.

A surviving bounded lift from [QFL2](BOUNDED_REALIZATION.md) must
therefore change the scalar order profile, avoid the scalar reduction,
or otherwise leave the stated hypotheses. Proving an order-profile
constraint for all possible bounded representatives would be the
additional all-degree theorem needed to turn this into a complete
nonlifting result. Searching increasingly large coefficient boxes
inside the (2,3) profile is unnecessary.

[Exact symbolic checks](check_curve_and_module.py) verify the general
Ore coefficient identities, the normalization, both indicial
factorizations and the two sharpness examples. The universal
algebraic-curve theorem uses the written Puiseux and
Riemann--Hurwitz proof, not a finite cover census. No formal
verification, external review, independent theorem replay, novelty
priority or DC2 solution is claimed.
