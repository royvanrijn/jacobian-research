# Oriented Davenport Cox maps

The Davenport boundary factorization admits a polynomial oriented lift which
removes the reciprocal derivative coordinate.  The result is a pair of
finite etale degree-seven Cox morphisms with constant residue Jacobian one,
one common oriented target, three source boundary components, and the full
point/line `GL_3(F_2)` Sunada monodromy.

This combines the Davenport--Sunada input with the multi-boundary Cox
ledger.  The resulting smooth charts are not affine spaces.

Work over

\[
 K=\mathbb Q(a),\qquad a^2+a+2=0.
\]

## 1. Boundary factorization

Use the Davenport polynomials `g_T(Y),h_T(Z)`, their common branch cubic
`Delta(T,u)`, and the factors from the
[boundary audit](DAVENPORT_COX_BOUNDARY_OBSTRUCTION.md):

\[
 \Delta(T,g_T(Y))=E_{3,g}E_{6,g}J_g^2,\qquad
 J_g=g_T'(Y),                                        \tag{1}
\]

\[
 \Delta(T,h_T(Z))=E_{3,h}E_{6,h}J_h^2,\qquad
 J_h=h_T'(Z).                                        \tag{2}
\]

The three factors on either source are distinct primes, with boundary
multiplicity vector `(1,1,2)`.

## 2. Oriented source and target

Define the common oriented target surface

\[
 \widetilde B=
 \{D^2=\Delta(T,u)\}\subset\mathbb A^3_{T,u,D}.       \tag{3}
\]

Define the two oriented source surfaces

\[
 \widetilde X_g=
 \{W^2=E_{3,g}(T,Y)E_{6,g}(T,Y)\}
 \subset\mathbb A^3_{T,Y,W},                         \tag{4}
\]

\[
 \widetilde X_h=
 \{W^2=E_{3,h}(T,Z)E_{6,h}(T,Z)\}
 \subset\mathbb A^3_{T,Z,W}.                         \tag{5}
\]

Restrict all three hypersurfaces to their smooth loci.  These opens contain
the generic points of each displayed boundary component.

Equations (1)--(2) define polynomial morphisms

\[
 \boxed{
 \widetilde\pi_g(T,Y,W)
 =(T,g_T(Y),D=J_gW),
 }                                                    \tag{6}
\]

\[
 \boxed{
 \widetilde\pi_h(T,Z,W)
 =(T,h_T(Z),D=J_hW).
 }                                                    \tag{7}
\]

Indeed,

\[
 (J_gW)^2=J_g^2E_{3,g}E_{6,g}
 =\Delta(T,g_T(Y)),
\]

and similarly for `h`.

## 3. Constant residue Jacobian

On `D!=0`, the target hypersurface residue form is

\[
 \Omega_{\widetilde B}
 =\frac{dT\wedge du}{2D}.                             \tag{8}
\]

On `W!=0`, use

\[
 \Omega_{\widetilde X_g}
 =\frac{dT\wedge dY}{2W},\qquad
 \Omega_{\widetilde X_h}
 =\frac{dT\wedge dZ}{2W}.                             \tag{9}
\]

Since

\[
 dT\wedge d(g_T(Y))=J_g\,dT\wedge dY,
\]

pullback through (6) gives

\[
 \widetilde\pi_g^*\Omega_{\widetilde B}
 =
 \frac{J_g\,dT\wedge dY}{2J_gW}
 =\Omega_{\widetilde X_g}.                           \tag{10}
\]

The same calculation applies to `h`.  Both sides extend as regular
generators across the generic oriented boundary, proving:

### Theorem 3.1

The maps (6)--(7) are finite etale morphisms of degree seven on their smooth
loci and have constant hypersurface-residue Jacobian one.

Finiteness follows because `Y` or `Z` satisfies a monic degree-seven
equation over the target and `W` is integral over that extension.  The
generic target orientation determines `W=D/J`, so adjoining the Cox square
root does not change the degree.

## 4. Three reduced boundary components

On `X_tilde_g`, the pullback of the oriented target parameter is

\[
 D=J_gW.
\]

The divisor `J_g=0` gives one component.  The divisor `W=0` has two prime
components, lying generically over `E_(3,g)=0` and `E_(6,g)=0`; locally the
source equations have the form

\[
 W^2=t\cdot(\text{unit}),
\]

so `W` is a uniformizer on the oriented source.

Consequently

\[
 \operatorname{div}(\widetilde\pi_g^*D)
 =
 \mathcal E_{3,g}+\mathcal E_{6,g}+\mathcal J_g,      \tag{11}
\]

with every coefficient equal to one.  The doubled derivative component in
the unoriented ledger `(1,1,2)` has been converted to the reduced oriented
ledger `(1,1,1)`.  The same statement holds for `h`.

Unlike the ordered-linear-factor charts, these covers are finite: the three
components are affine boundary divisors over `D=0`, not dicritical divisors
at source infinity.  The orientation solves the determinant problem but
does not create nonproperness.

## 5. Sunada monodromy

Away from `D=0`, the oriented construction is a finite base change of the
original Davenport point and line covers.  Their common Galois closure
therefore remains the `GL_3(F_2)` cover, and the two degree-seven actions
remain the nonconjugate point and line actions with identical permutation
characters.

Thus (6)--(7) give a pair of nonisomorphic constant-Jacobian Cox covers over
one common target, with identical zeta functions on every common good
oriented fiber.

The group contains elements with zero, one, three, and seven fixed points in
the degree-seven action.  Hence this global cover is not a permutation
family over all rational target points.  The specialization `T=0` degenerates
to the toric seventh-power map and can be a permutation when
`gcd(7,q-1)=1`, but that specialization no longer carries the generic
`GL_3(F_2)` monodromy.

This cannot be repaired by an arithmetic twist preserving the geometric
degree-seven action.  The
[exceptional-twist obstruction](DAVENPORT_EXCEPTIONAL_TWIST_OBSTRUCTION.md)
proves

\[
N_{S_7}(GL_3(\mathbb F_2))=GL_3(\mathbb F_2),
\]

so arithmetic monodromy cannot be a larger exceptional extension of the
same point or line action.

## 6. What is now complete and what is not

For the Davenport pair, the following pieces are complete:

1. the full height-one pullback `(1,1,2)`;
2. the source and target unit lattices;
3. the reciprocal one-coordinate Cox suspension;
4. the polynomial oriented lift `(1,1,1)`;
5. constant residue Jacobian one;
6. one common target and one common Galois closure; and
7. exact point/line fiber-zeta equality.

The remaining absolute problem is geometric, not group-theoretic: fill or
modify the three source boundary components and the target boundary so that
both sides become affine spaces without destroying (10) or the common
monodromy.  Stable straightening and coordinate-preserving polynomial
suspensions are already ruled out by the boundary audit.

## 7. Reproduction

Run

```bash
.venv/bin/python scripts/verify_oriented_davenport_cox_maps.py
```

The checker verifies both boundary factorizations, the oriented square
relations, the residue cancellation, the reduced three-component ledger,
the degree-seven point/line actions, and the fixed-point obstruction to a
global permutation family.
