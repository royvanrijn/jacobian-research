# Positive invariance regression for the foundational Keller map

Most boundary checks compare examples expected to differ.  This audit tests
the opposite requirement: hostile presentations of the same foundational Keller map must
produce canonically isomorphic boundary objects.

The regression simultaneously changes source coordinates, target
coordinates, primitive element, projective root coordinate, and dimension.
Raw equations are intentionally different.  The comparison is made through
the transported finite function-field extension and normalization.

## 1. Original presentation

Write the foundational Keller map as `F=(a,b,c)` and put

\[
 P(T)=cT^3-2T^2+bT-2a.                                    \tag{1}
\]

The marked root is `t=y+1/x` on the first reconstruction chart.  The
canonical normalization is the normalization of the finite marked-root
incidence, and its boundary consists of one prime over the irreducible target
discriminant

\[
 \Delta(a,b,c)=\operatorname{Disc}_T(P).                   \tag{2}
\]

The boundary prime maps birationally with `(e,f)=(2,1)`.  Thus every level of
the invariant ladder has one labelled target vertex and one upstairs vertex;
the formal level additionally retains their completed finite map and the
different.

## 2. Hostile source and target coordinates

Use new source coordinates `(X,Y,Z)` and the triangular automorphism

\[
 \alpha(X,Y,Z)=
 \left(X+Y^2+Z,\ Y+Z^2,\ Z\right).                         \tag{3}
\]

Use the target automorphism

\[
 \beta(a,b,c)=
 \left(a+b^2+c,\ b+c^2,\ c\right),                         \tag{4}
\]

whose inverse is

\[
 c=C,\qquad b=B-C^2,\qquad
 a=A-(B-C^2)^2-C.                                          \tag{5}
\]

Set

\[
 G=\beta\circ F\circ\alpha.                               \tag{6}
\]

Both triangular automorphisms have Jacobian one, so `det DG=-2`.  More
importantly, (3)--(5) identify the function-field inclusions for `F` and `G`.
Functoriality of integral closure therefore gives an isomorphism of canonical
triples

\[
 (X\subset\overline X_F\to Y)
 \simeq
 (X'\subset\overline X_G\to Y').                           \tag{7}
\]

The transformed target boundary is not compared by its displayed equation:
it is exactly `beta(V(Delta))`, with equation obtained by substituting (5)
in (2).  The executable independently verifies that this is also the
discriminant of the transformed cubic resolvent.

## 3. A different primitive element

Replace the marked root `T` by

\[
 \Theta=T+T^2.                                             \tag{8}
\]

Eliminating `T` produces the visibly different cubic

\[
 M(\Theta)=\operatorname{Res}_T(P(T),\Theta-T-T^2).         \tag{9}
\]

This primitive-element change does not change the field.  Indeed, modulo
`P(T)` one obtains

\[
 T=
 { (c+2)\Theta+2a
  \over c(\Theta+1)+b+2}.                                  \tag{10}
\]

Thus `K(T)=K(Theta)` on the generic fiber.  The normalizations of (1) and
(9) are the same integral closure.  In fact the audit finds the deliberately
different raw discriminant

\[
 \operatorname{Disc}_\Theta(M)
 =\operatorname{Disc}_T(P)
 \{2ac-bc-2b-c^2-4c-4\}^2.
\]

The extra square is a primitive-generator collision locus, not a new branch
divisor of the normalized field extension.  This is precisely the kind of
coordinate dependence the positive regression is meant to catch.  The test
compares (9)--(10) and the integral closures, not the raw discriminants.

## 4. A different compactification coordinate

On the marked projective root line use

\[
 V={T+1\over T+2},qquad T={2V-1\over1-V}.                 \tag{11}
\]

The transformed binary cubic is

\[
 \widetilde P(V)=(1-V)^3
 P\left({2V-1\over1-V}\right).                             \tag{12}
\]

The Möbius matrix in (11) has determinant one.  Hence (12) defines an
isomorphic finite projective incidence, carries the repeated-root divisor to
the repeated-root divisor, and has the same discriminant (up to the expected
unit for a general matrix).  This explicitly tests that changing root charts
or compactification coordinates does not change the canonical boundary.

## 5. Stabilization

Adjoin an identity coordinate `W` and compare

\[
 F\times\operatorname{id}_{\mathbb A^1},qquad
 G\times\operatorname{id}_{\mathbb A^1}.                   \tag{13}
\]

The finite normalization, boundary prime, discriminant divisor, reduced and
scheme intersections, completed maps, different, conductor data, and
valuation filtration are their former counterparts times `A^1` or the
corresponding completed cylinders.  The checker verifies the block Jacobian
and that the transported discriminant ideal is extended from the
three-dimensional target.

## 6. Canonical output comparison

The regression compares the following maps, not strings of equations:

\[
 \begin{array}{c}
 k(a,b,c)\subset k(x,y,z)\\
 \downarrow\scriptstyle\beta^*\qquad\downarrow\scriptstyle\alpha^*\\
 k(A,B,C)\subset k(X,Y,Z),
 \end{array}                                               \tag{14}
\]

together with the generic-field identifications

\[
 K(T)=K(\Theta)=K(V).                                      \tag{15}
\]

Integral closure, the distinguished affine open, boundary primes,
intersections, completions, differents, finite-stratum conductors, and
valuation filtrations are functorial under these maps.  Consequently

\[
 \mathfrak I^*(F\times\operatorname{id})
 \simeq
 \mathfrak I^*(G\times\operatorname{id}),
 \qquad *=\mathrm{red},\mathrm{sch},\mathrm{formal}.       \tag{16}
\]

## Independent executable

Run

```bash
.venv/bin/python scripts/audit_foundational_invariance_regression.py
```

It verifies exactly:

- invertibility and Jacobians of (3)--(6);
- the transported marked-root identity;
- transport of the discriminant prime under (4)--(5);
- the resultant resolvent (9) and rational inverse (10);
- the Möbius chart identity and discriminant covariance; and
- the stabilized block Jacobian and extended boundary ideal.

This is a positive covariance test: every presentation must return the same
canonical object even though none of the displayed resolvents or coordinate
equations is required to match literally.
