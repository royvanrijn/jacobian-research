# Exact geometric-degree spectrum in dimension three

For a dominant polynomial map

\[
F:\mathbb A^3_{\mathbb C}\longrightarrow\mathbb A^3_{\mathbb C},
\]

write

\[
\operatorname{gdeg}(F)
=\bigl[\mathbb C(x_1,x_2,x_3):
         \mathbb C(F_1,F_2,F_3)\bigr].
\]

For a Keller map this is also the number of points in a generic geometric
fiber: the nonzero constant Jacobian makes the induced function-field
extension finite and separable.

## Spectrum theorem

The exact set of geometric degrees of noninvertible Keller maps
`A^3_C -> A^3_C` is

\[
\boxed{\{3,4,5,\ldots\}}.
\]

### Exclusion of degrees one and two

If `gdeg(F)=1`, then the induced function-field extension is trivial, so `F`
has a rational inverse. The classical birational case of the Jacobian theorem
makes this inverse polynomial; hence `F` is a polynomial automorphism.

If `gdeg(F)=2`, then

\[
\mathbb C(x_1,x_2,x_3)/\mathbb C(F_1,F_2,F_3)
\]

is a separable quadratic extension and is therefore Galois. Campbell's
Galois theorem makes `F` a polynomial automorphism. Razar and independently
Wright gave subsequent algebraic treatments.

The same exclusion works over every characteristic-zero ground field. The
[finite-étale Keller-fiber note](FINITE_ETALE_KELLER_FIBERS.md) gives the
explicit descent: move the coefficients to a finitely generated subfield,
embed it into `C`, use stability of generic degree under scalar extension,
and descend the unique formal/polynomial inverse. Thus no
algebraic-closedness assumption is required.

Every noninvertible Keller map therefore has geometric degree at least three.

### Realization of every degree at least three

The foundational map realizes degree three. Uniformly, for every `N>=3` the
split polynomial

\[
H_N(W)=W^2(1-W)(1+W^{N-3})
\]

has exact degree `N` and satisfies the weighted admissibility conditions

\[
H_N(0)=H_N'(0)=H_N(1)=0,\qquad
H_N'(1)=-2,\qquad
\frac{H_N''(1)}2=-(N+1)\ne-2.
\]

Taking `c=2` and `b_0=1` in the weighted construction therefore gives a
polynomial map `G_{H_N}:A^3_C -> A^3_C` with

\[
\det DG_{H_N}=2.
\]

Its generic inverse equation is

\[
H_N(W)-BCW+2AC^2=0.
\]

On the generic locus `C!=0`, every simple root reconstructs one source point,
so the tangent-map and weighted-suspension theorem gives
`gdeg(G_{H_N})=N`. Since `N>1`, this map cannot be invertible. This realizes
every integer in the displayed spectrum. At `N=3`, `H_3=2W^2(1-W)` is the
foundational seed up to the harmless normalization used in the foundational
model.

There is now a stronger arithmetic realization route as well. The
[finite-étale Keller-fiber theorem](FINITE_ETALE_KELLER_FIBERS.md) constructs,
from every squarefree polynomial of degree `N>=3`, a determinant-minus-two
map of `A^3` having its finite étale algebra as one full fiber. In particular,
it gives another direct realization of every degree in the spectrum.

## Left--right invariance

If

\[
F'=L\circ F\circ R
\]

for polynomial automorphisms `L` and `R`, their induced function-field
isomorphisms identify the two extensions defining geometric degree. Hence

\[
\operatorname{gdeg}(F')=\operatorname{gdeg}(F).
\]

Adjoining identity coordinates only adds the same algebraically independent
variables to the top and bottom fields, so geometric degree is also invariant
under stable polynomial left--right equivalence. In particular the maps
realizing different `N` above are automatically left--right, and even stably
left--right, inequivalent.

## Classical input

- L. Andrew Campbell,
  [*A condition for a polynomial map to be invertible*](https://doi.org/10.1007/BF01349234),
  Math. Ann. 205 (1973), 243--248.
- Michael Razar,
  [*Polynomial maps with constant Jacobian*](https://doi.org/10.1007/BF02764906),
  Israel J. Math. 32 (1979), 97--106.
- David Wright,
  [*On the Jacobian conjecture*](https://doi.org/10.1215/ijm/1256047158),
  Illinois J. Math. 25 (1981), 423--440.

The internal existence and degree inputs are the
[foundational geometry](FOUNDATIONAL_GEOMETRY.md), the
[tangent-map core](TANGENT_MAP_CORE.md), the
[weighted-seed theorem](WEIGHTED_SEED_THEOREM.md), and the
[finite-étale Keller-fiber theorem](FINITE_ETALE_KELLER_FIBERS.md).
