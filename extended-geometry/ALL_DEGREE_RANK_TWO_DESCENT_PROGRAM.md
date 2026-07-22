# Toward rank-two symplectic descent in every degree

This note defines a workable obstruction, proves rank-two completion on the
normalized degree-five seed surface, and gives the next all-degree
computations.  The degree-five result is an internal repository theorem with
an exact two-parameter certificate and no recorded external specialist review.

## 1. Normalized seeds and fixed-third-component slices

A normalized degree-`N` weighted seed is

\[
 H(W)=\sum_{j=3}^N h_j(W^j-W^2),\qquad H'(1)=-1.
\]

It has dimension `N-3`.  Put

\[
 \kappa=H''(1),\qquad
 a=-\frac{1+\kappa}{2+\kappa}.
\]

The third component of the weighted map is

\[
 C=x(1+a,xy+x^2z).
\]

Thus fixing `gamma`, after fixing the harmless nonzero coefficient of `x^2z`,
is exactly the linear equation `H''(1)=kappa`.  It cuts the normalized seed
space from dimension `N-3` to `N-4`, provided this functional is nonconstant.
This proves the proposed dimension count for the *seed slice*; it does not by
itself prove completion or inequivalence of the resulting four-dimensional
maps.

For `a != 0`, every such slice has the same adapted quotient coordinate.  Set

\[
 W=Z+\psi(Q),\qquad Y=Q-\frac{XW}{3},
\]

and make the seed-dependent diagonal source change

\[
 (x,y,z)=\left(X,-\frac{3}{2a}Y,-\frac12W\right).
\]

Then, identically and independently of the seed and of `psi`,

\[
 2C=R=2X-3X^2Q.
\]

This is the common rank-two base needed in every degree.  The value `a=0`
(`kappa=-1`) requires another adapted chart; `kappa=-2` is outside the
weighted construction.

## 2. A rank-two flux obstruction

Let `K` be the rational function field of a bounded seed family and put

\[
 B=K[X,Q,Z],\qquad B_X=B[X^{-1}].
\]

On `B` use the quotient Poisson bracket

\[
 \{X,Z\}=-3X^2,\qquad \{Q,Z\}=6XQ-2,\qquad \{X,Q\}=0.
\]

For an adapted triple `Phi=(S,T,R)` with

\[
 \{S,T\}=1,\qquad \{R,S\}=\{R,T\}=0,
\]

let `w_Phi` be the horizontal derivation characterized by

\[
 w_\Phi(S)=w_\Phi(T)=0,\qquad w_\Phi(R)=1.
\]

Let `w_E={E,-}` be the fixed derivation in the common four-dimensional
coordinate system.  The difference

\[
 \delta_\Phi=w_\Phi-w_E
\]

has a normalized Hamiltonian `f_Phi in B_X`: integrate its `X`-component in
`Z` with zero `Z`-constant, and then integrate the remaining component after

\[
 v=X^{-1},\qquad \rho=R,qquad Q=\frac{v(2-\rho v)}3.
\]

This is the finite homotopy already used in the degree-five proof.  Define the
raw residue

\[
 \operatorname{Res}_X(\Phi)=\operatorname{PP}_{X=0}(f_\Phi),
\]

the sum of the negative powers of `X`.  For seeds and shears of bounded
degree, denominator and numerator degrees are uniformly bounded, so these
residues lie in a finite-dimensional Laurent space `V_{N,d}`.

Let `Gamma_{R,d}` be the bounded-degree group of base automorphisms preserving
`R`; it contains `Z -> Z+psi(Q)`.  The proposed obstruction is the orbit

\[
 \boxed{\operatorname{Flux}_{N,d}(H)
 = [\operatorname{Res}_X(\Phi_H)]\in V_{N,d}/\Gamma_{R,d}.}
\]

This quotient should initially be treated as an orbit set (or quotient
stack), not automatically as a vector-space quotient.  On an ansatz where
the shear response is linear, it can be replaced by the cokernel of the
linear response map.  The zero-orbit criterion is exact:

\[
 \operatorname{Flux}_{N,d}(H)=0
 \quad\Longleftrightarrow\quad
 \text{an allowed base change makes }f_\Phi\text{ polynomial}.
\]

Constants and linear terms of `psi` belong to the elementary base gauge, so
the first genuine coefficient is the coefficient of `Q^2`.

## 3. The full normalized degree-five surface

The normalized degree-five surface has coordinates `(kappa,tau)`:

\[
 H_{\kappa,\tau}(W)=W^2(W-1)
 \left(
 \tau W^2+\left(\frac\kappa2-2\tau+2\right)W
 -\frac\kappa2+\tau-3
 \right).
\]

It satisfies `H'(1)=-1` and `H''(1)=kappa`.  It is useful in the exact
calculation to replace `kappa` by

\[
 a=-\frac{1+\kappa}{2+\kappa},\qquad
 \kappa=-\frac{1+2a}{1+a}.
\]

With `psi(Q)=s_2Q^2`, direct calculation over `Q(a,tau,s_2)` proves that the
localized Hamiltonian has denominator

\[
 15482880a^8(a+1)^2X^4
\]

and that its four low numerator coefficients are respectively

\[
 4096a^8\mathcal O,\quad
 24576Qa^8\mathcal O,\quad
 92160Q^2a^8\mathcal O,\quad
 276480Q^3a^8\mathcal O,
\]

where

\[
\begin{aligned}
\mathcal O={}&28(a+1)^2s_2+12a(a+1)^2\tau^2\\
&-18a(a+1)(4a+5)\tau
-216a^3-648a^2-738a-315.
\end{aligned}
\]

There are no lower powers of `X`.  Solving `mathcal O=0` and returning to
`kappa` gives

\[
 \boxed{
 s_2=K(\kappa,\tau):=
 \frac{
 12(\kappa+1)\tau^2
 -18(\kappa+1)(\kappa+6)\tau
 +9(\kappa^3+16\kappa^2+52\kappa+72)
 }{28(\kappa+2)}.}
\]

Up to a nonzero parameter unit, the raw residue is

\[
 (s_2-K(\kappa,\tau))
 \left(
 X^{-4}+6QX^{-3}+\frac{45}{2}Q^2X^{-2}
 +\frac{135}{2}Q^3X^{-1}
 \right).
\]

### Theorem — degree-five surface completion

On the chart `kappa != -2,-1`, the normalized Hamiltonian is polynomial
exactly for `s_2=K(kappa,tau)`.  The resulting four coordinates have the six
canonical brackets, and the source change is a polynomial automorphism.
Therefore every admissible weighted degree-five map on this chart has an
exact symplectic completion in `A^4`, polynomially left--right equivalent to
the original map times `id_(A^1)`.

Consequently:

* the completable locus on `kappa != -2,-1` is the entire two-dimensional
  degree-five seed surface, not a divisor;
* the completing quadratic shear is unique modulo constant and linear base
  gauge; and
* allowing higher-degree `psi` is unnecessary in degree five.

This correctly recovers the published line.  There `kappa=-9` and

\[
 \tau=\frac{3}{2(\lambda-1)},
\]

so the completing shear becomes

\[
 K\left(-9,\frac{3}{2(\lambda-1)}\right)
 =-\frac{27(57\lambda^2-138\lambda+73)}
 {196(\lambda-1)^2}.
\]

The published theorem is therefore the generic fixed-`gamma` line inside this
surface.  On the ordinary boundary-clean open where the decorated
normalization invariant is generically finite, completion transfers the full
two-dimensional stable moduli to exact symplectic maps of `A^4`.  Thus degree
five is stronger than the fixed-`gamma` target `N-4=1`: it gives `N-3=2`
dimensions on this chart.

The exact certificate is
[`verify_degree_five_flux_surface.py`](../scripts/verify_degree_five_flux_surface.py).
The earlier slice-oriented exploratory replay remains
[`explore_degree_five_flux_surface.py`](../scripts/explore_degree_five_flux_surface.py).

### Exceptional divisors and the moduli open

The exact calculation separates four kinds of exceptional parameters:

1. `tau=0` is the degree-drop divisor: the coefficient of `W^5` in `H` is
   `tau`.
2. `kappa=-2` is genuinely outside the weighted marked-root construction,
   because `a=-(1+kappa)/(2+kappa)` is undefined.
3. `kappa=-1` gives `a=0`.  The shear itself remains regular,
   `K(-1,tau)=45/4`, but the adapted source change divides by `a`.  This is a
   chart exception; the calculation does not prove non-completability there.
4. The ordinary-cusp condition removes the divisor

\[
\begin{aligned}
\Delta_{\rm cusp}(\kappa,\tau)={}&
\kappa^4-4\kappa^3\tau+18\kappa^3
+8\kappa^2\tau^2-36\kappa^2\tau+123\kappa^2\\
&-8\kappa\tau^3+24\kappa\tau^2-72\kappa\tau+376\kappa
+8\tau^4+24\tau^3-24\tau^2+72\tau+432.
\end{aligned}
\]

Indeed, `disc(H'')=432 Delta_cusp`.  The decorated-normalization theorem
removes additional proper incidence divisors for cusp-image collisions and
nonordinary bitangencies.  Their complement is already known to be a nonempty
open.  Intersecting it with

\[
 \tau(\kappa+1)(\kappa+2)\Delta_{\rm cusp}\ne0
\]

is therefore a nonempty two-dimensional open on which the maps have degree
five, the rank-two completion is defined, and the stable invariant is
generically finite.  No resolution of the chart divisor `kappa=-1` is needed
for the two-dimensional theorem.

## 4. Degree-by-degree computation

For each `N>=5`:

1. Parametrize normalized seeds by `N-3` affine coordinates and use `kappa`
   as one coordinate.
2. Make the common diagonal change above on `a != 0` and compute the adapted
   triple `(S,T,R)` without expanding the source map first.
3. Choose `psi(Q)=sum_{j=2}^d s_jQ^j`, initially with `d=N-3`, and construct
   the normalized localized Hamiltonian by the finite homotopy.
4. Extract every coefficient of `X^{-1},...,X^{-M}`.  Clear only parameter
   denominators, not powers of `X`.
5. Eliminate the shear coefficients `s_j`.  The elimination ideal in seed
   parameters defines the completable locus; its graph records all completing
   shears.
6. Repeat after imposing `kappa=kappa_0`.  Compare the rank of the shear
   response and the dimension of the resulting locus with `N-4`.
7. Certify stable moduli only after restricting the decorated-normalization
   invariant to the completable locus and proving that its generic fiber is
   finite there.

The degree-five result warns against assuming that fixed `gamma` contributes
an additional completion equation: at least in degree five, `kappa` varies and
the shear absorbs the flux throughout the surface.

## 5. Quantization to the second Weyl algebra

Once a polynomial exact symplectic map

\[
 G=(R,T,D,S):\mathbb A^4\to\mathbb A^4
\]

is available, a direct lift to the second Weyl algebra must solve more than
the associated-graded problem.  Choose Weyl generators with

\[
 [D,R]=[S,T]=\hbar
\]

and all mixed commutators zero.  A filtered quantization seeks ordered lifts
`R_hbar,T_hbar,D_hbar,S_hbar` with the exact same relations.  The principal
symbols give `G`, but at filtration order `m-2` the commutator differs from
the Poisson bracket by an ordering cocycle.  Corrections must therefore be
found recursively from linear equations of the form

\[
 d_G(\Delta_m)=-\mathcal O_m,
\]

where `d_G` is the linearized commutator map and `O_m` is determined by lower
orders.  Nonvanishing classes in the corresponding cokernel are genuine
ordering obstructions.

The practical first test is the completed degree-five family: use symmetric
(Weyl) ordering, compute the six commutator defects through the first two
lower filtration levels, and solve simultaneously for corrections.  Success
would give a filtered endomorphism of `A_2`; injectivity follows from the
injective associated-graded map, while non-surjectivity would still require a
separate filtered argument.
