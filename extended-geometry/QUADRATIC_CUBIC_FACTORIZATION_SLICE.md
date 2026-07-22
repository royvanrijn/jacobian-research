# The normalized quadratic-cubic factorization slice

This note studies the first consecutive-degree case beyond the foundational
linear-quadratic slice.  The conclusion is negative but informative: the
natural normalized source is a smooth factorial affine fivefold with trivial
canonical class, yet it is not affine five-space.

Work over a field `k` of characteristic zero.  Write

\[
 A=a_0T^2+a_1TS+a_2S^2,\qquad
 B=b_0T^3+b_1T^2S+b_2TS^2+b_3S^3
\]

and use the tangent coefficient

\[
 m=[AB]_{T^4S}=a_0b_1+a_1b_0.
\]

The normalized affine slice is

\[
 X_{2,3}=\{m=1,\ \operatorname{Res}(A,B)=1\}\subset\mathbb A^7. \tag{1}
\]

The unequal-degree theorem makes multiplication etale on (1).  Consecutive
degrees make the normalization unique, so projectivization identifies (1)
with

\[
 U_{2,3}=(\mathbb P^2\times\mathbb P^3)\setminus(R\cup E), \tag{2}
\]

where `R` is the resultant divisor and `E={m=0}`.  A squarefree binary
quintic has `binomial(5,2)=10` quadratic factors, so multiplication is
generically ten-to-one.

All conclusions below concern this natural tangent coefficient.  Other
hyperplanes in the dual binary-quintic space need not have the same
intersection geometry.

## 1. Units and Picard group

The ambient Picard group is

\[
 \operatorname{Pic}(\mathbb P^2\times\mathbb P^3)=\mathbb Z^2.
\]

The binary resultant is irreducible of bidegree `(3,2)`, while the rank-two
bilinear form `m` is irreducible of bidegree `(1,1)`.  Divisor localization
therefore uses the matrix

\[
 \begin{pmatrix}3&1\\2&1\end{pmatrix},\qquad \det=1.  \tag{3}
\]

It is an isomorphism of lattices.  Hence

\[
 \boxed{\mathcal O(U_{2,3})^*/k^*=0,\qquad
 \operatorname{Pic}(U_{2,3})=0.}                      \tag{4}
\]

Thus the first obstructions agree with `A^5`; they do not decide the
question.

## 2. Grothendieck class

Let `L=[A^1]` and `[P^j]=1+L+...+L^j`.  Stratify the resultant divisor by the
degree of `gcd(A,B)`.

- For gcd degree one, choose the common linear form and a coprime residual
  `(1,2)` pair.  The latter has class
  \[
   [P^1][P^2]-[P^1]^2=L^2+L^3.
  \]
- For gcd degree two, choose the quadratic gcd and the remaining linear
  factor, giving `P^2 x P^1`.

Consequently

\[
 [\{\operatorname{Res}\ne0\}]=L^4+L^5.               \tag{5}
\]

It remains to impose `m!=0`.  The coprime part of `E` has class `L^4+L^3`.
Here is a direct stratification.

If `a_0!=0`, make `A` monic.  Polynomial division identifies the hyperplane
of possible `B` with triples `[b_0:r_1:r_0]`, where
`r_1T+r_0` is the remainder modulo `A`.  Coprimality says this nonzero linear
remainder does not vanish at a root of `A`.  Projecting to its root gives a
space of class `L^3`, and the free `b_0` contributes another `L`; this chart
has class `L^4`.

If `a_0=0`, the equation `m=0` and coprimality force `A=S^2` and `b_0!=0`.
The remaining cubic coefficients form `A^3`, contributing `L^3`.  Therefore

\[
 [E\cap\{\operatorname{Res}\ne0\}]=L^4+L^3.          \tag{6}
\]

Subtracting (6) from (5) gives

\[
 \boxed{[U_{2,3}]=L^5-L^3=L^2[\mathrm{SL}_2].}        \tag{7}
\]

The Hodge--Deligne realization sends this to

\[
 E_c(U_{2,3};u,v)=(uv)^5-(uv)^3,\qquad
 E(U_{2,3};u,v)=1-(uv)^2.                              \tag{8}
\]

In particular (7)--(8) prove

\[
 \boxed{U_{2,3}\not\simeq\mathbb A^5}.               \tag{9}
\]

The suggestive equality with `L^2[SL_2]` is motivic only; no product
isomorphism `U_{2,3}=A^2 x SL_2` is asserted.

## 3. Finite-field counts

The same strata are defined over every finite field and give

\[
 \boxed{\#U_{2,3}(\mathbb F_q)=q^5-q^3.}              \tag{10}
\]

This differs from `#A^5(F_q)=q^5` for every prime power.  It is an independent
realization of the motivic obstruction and avoids relying on cancellation in
the Grothendieck ring.

## 4. Low-degree cohomology: what is determined

Over `C`, the standard meridian sequence for the complement of two prime
divisors begins

\[
 0=H^1(\mathbb P^2\mathbin\times\mathbb P^3,\mathbb Q)
 \longrightarrow H^1(U_{2,3},\mathbb Q)
 \longrightarrow \mathbb Q(-1)^2
 \xrightarrow{\delta}H^2(\mathbb P^2\mathbin\times\mathbb P^3,\mathbb Q).
\]

In the two hyperplane bases, `delta` is exactly the divisor-class matrix
(3).  It is invertible.  The open is connected, so this proves

\[
 \boxed{H^0(U_{2,3},\mathbb Q)=\mathbb Q,\qquad
 H^1(U_{2,3},\mathbb Q)=0.}                          \tag{11}
\]

Equation (8) computes the complete **virtual** mixed Hodge polynomial.  Its
first nonconstant term is `-(uv)^2`, as for `A^2 x SL_2`.  It is tempting to
identify this with `H^3=Q(-2)` and all other positive-degree rational
cohomology zero.  That conclusion does not follow from the Grothendieck class:
different cohomological degrees can cancel in the Hodge--Deligne polynomial.
A full normal-crossings or simplicial resolution of `R union E`, including
the degree-two gcd conductor, is still required to separate `H^2` and `H^3`.
Accordingly the proved low-degree statement here is (11), not the stronger
`SL_2`-shaped cohomology guess.

For context, Farb--Wolfson compute the topology and arithmetic of the
`resultant=1` space for **two monic forms of the same degree** in
[`Topology and arithmetic of resultants, II`](https://arxiv.org/abs/1507.01283).
Their purity theorem is useful guidance, but it does not directly cover this
unequal-degree projective complement with the additional tangent divisor.
No cohomology assertion above is imported from that paper.

## 5. Factoriality and canonical class

The slice is smooth because it is the base change of the etale
coefficient-resultant map.  Hence it is regular and locally factorial.  For a
regular affine variety, `Cl=Pic`; (4) therefore proves

\[
 \boxed{k[X_{2,3}]\text{ is a UFD}.}                  \tag{12}
\]

The ambient canonical class in the projective model is `(-3,-4)`.  In the
boundary lattice,

\[
 (-3,-4)=(3,2)-6(1,1).                                \tag{13}
\]

Both boundary classes vanish after restriction to `U`, so

\[
 \boxed{K_{U_{2,3}}=0.}                               \tag{14}
\]

Equivalently, (1) is a smooth codimension-two complete intersection in
affine seven-space, and the Poincaré-residue form trivializes its canonical
bundle.

## 6. Consequence

The natural `(2,3)` multiplication map is an etale, generically degree-ten
map

\[
 U_{2,3}\longrightarrow\mathbb A^5,
\]

but its source is not `A^5`.  It therefore does **not** produce a new `JC(5)`
counterexample.  It is still a genuinely new nonproper etale cover model and
its `SL_2`-shaped virtual motive suggests two next questions: compute the full
boundary-resolution cohomology, and test whether a different binary-quintic
hyperplane can remove the motivic obstruction.

The exact class calculation and direct finite-field enumeration are checked
by
[`verify_quadratic_cubic_factorization_invariants.py`](../scripts/verify_quadratic_cubic_factorization_invariants.py).
