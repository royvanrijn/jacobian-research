# Fixed-rank Dixmier reduction: the first invariant-pair search

The repository has two different lower-rank constructions which must not be
conflated.

1. The inverse-Jacobian construction quantizes every three-variable Keller
   map as an injective non-surjective endomorphism of the third Weyl algebra.
2. The rank-two residue calculation gives a noninjective polynomial Poisson
   endomorphism on four variables, with two canonical pairs, which is
   polynomially left--right equivalent to the original Keller map times one
   identity coordinate.

The second statement is already a classical rank-two compression.  It is not
an endomorphism of the second Weyl algebra.  Its direct global filtered
quantization is precisely the missing step: the standard parity-preserving
ansatz is obstructed at `hbar^5`, while the formal-local quantization does not
control polynomiality or descent across `X=0`.

This note starts a complementary search inside the explicit `A_3`
endomorphism itself.

There is also a first-class localized route.  The foundational map is
equivariant for a natural `G_m` action.  Its full affine moment reductions
are singular, but the primitive target-momentum chart has target quotient
`A^4`, a smooth explicit source hypersurface, and retains all three collision
points as distinct quotient points.  The remaining questions are polynomial
triviality of that source fourfold and its localized quantum normalizer
quotient.  See
[Natural torus Hamiltonian reduction](NATURAL_TORUS_HAMILTONIAN_REDUCTION.md).

## 1. What an invariant pair would have to prove

Write the foundational Keller map as `F` and its inverse-Jacobian Weyl lift as

\[
 \Phi_F(X_i)=F_i(X),\qquad
 \Phi_F(D_i)=\sum_r(DF^{-1})_{ri}(X)D_r.
\]

A useful spectator pair for reduction would be a copy

\[
 B=k\langle h,L\rangle\simeq A_1,
 \qquad [L,h]=1,
\]

which is preserved by `Phi_F`.  For descent to `A_2`, invariance by itself is
not sufficient.  One needs one of the following stronger structures.

* **Centralizer/tensor-factor reduction:** prove
  `A_3 = B tensor C` (or a suitable double-centralizer statement) with
  `C` isomorphic to `A_2`, and prove `Phi_F(C)` is contained in `C`.
* **Quantum Hamiltonian reduction:** choose a moment element `mu` for the
  pair, prove `Phi_F` preserves the normalizer of the constraint
  `A_3(mu-c)`, and identify the normalizer quotient with `A_2`.
* **Filtered corner or module reduction:** construct a faithful reduced
  algebra and prove that its associated graded is the desired rank-four
  symplectic leaf.

The Weyl algebra is simple, so there is no nonzero two-sided ideal obtained by
merely setting one Weyl generator equal to a scalar.  A quotient must use a
Hamiltonian normalizer, a module/corner, or a centralizer construction.

There is also a direction issue.  If the restriction `Phi_F|B` were itself
injective and non-surjective, it would already be a rank-one Dixmier
counterexample.  For a rank-two descent, the more natural spectator behavior
is an automorphic (ideally trivial) action on `B`, with non-surjectivity
proved on the complementary reduced algebra.

## 2. Cotangent graphs and the intrinsic graph--flux criterion

The four-dimensional classical completion has a universal part which is
useful independently of its explicit formulas.  Let

\[
 F_0=(S,T,R):\mathbb A^3\longrightarrow\mathbb A^3
\]

be any Keller triple.  On `A^3 times A^1_lambda` put

\[
 \Omega_F=dS\wedge dT+d\lambda\wedge dR.             \tag{GF1}
\]

This is a polynomial symplectic form: its square is a nonzero constant
multiple of `dS wedge dT wedge d(lambda) wedge dR`.  The cotangent-graph
embeddings from the Poisson audit identify

\[
 (F_0\times\operatorname{id})^*(ds\wedge dt+d\lambda\wedge dr)
 =\Omega_F.                                           \tag{GF2}
\]

Consequently the exact rank-drop problem is the following.

> **Graph--Darboux criterion.**  The ordered triple `(S,T,R)` has a
> polynomial rank-two completion if and only if there is a polynomial
> automorphism
> \[
> H:\mathbb A^4_{\rm std}\longrightarrow
>   \mathbb A^3\times\mathbb A^1_\lambda
> \]
> with `H^* Omega_F=omega_std`.  In that case
> \[
> G=(R,T,\lambda,S)\circ H                           \tag{GF3}
> \]
> is a polynomial symplectic endomorphism, and it has the same fiber schemes
> as `F_0` after adjoining one identity coordinate.

This criterion is coordinate-free once the ordered target polarization
`(S,T,R)` is fixed.  The obstruction is not algebraic de Rham cohomology:
`Omega_F` is exact.  It is the polynomial-automorphism orbit

\[
 \boxed{
 \operatorname{GD}(F;S,T,R)
 =[\Omega_F]\in
 \{\text{polynomial symplectic forms on }\mathbb A^4\}/
 \operatorname{Aut}_{\rm poly}(\mathbb A^4).}         \tag{GF4}
\]

The graph--Darboux obstruction vanishes when this orbit contains the
standard form.

On an adapted chart with quotient algebra `B=k[X,Q,Z]`, Casimir `R`, and
slice `E` satisfying `{E,R}=1`, the orbit problem becomes computable.  Let
`w_(S,T,R)` be the horizontal derivation killing `S,T` and sending `R` to
one, and let `w_E={E,-}`.  After localizing along the chosen boundary
equation `X`, solve

\[
 w_{(S,T,R)}-w_E=\{f,-\},\qquad f\in B[X^{-1}].       \tag{GF5}
\]

Normalize `f` by the finite homotopy used in the rank-two descent.  Its
negative `X`-principal part, modulo the action of bounded `R`-preserving
Poisson base automorphisms, is the **graph--flux class**

\[
 \operatorname{Flux}_X(F;S,T,R)
 =[\operatorname{PP}_{X=0}(f)].                      \tag{GF6}
\]

For the declared automorphism class this is an exact obstruction:

\[
 \operatorname{Flux}_X=0
 \Longleftrightarrow
 \text{an allowed base change makes }f\text{ polynomial}
 \Longleftrightarrow
 D=E+f\text{ completes the Darboux system}.          \tag{GF7}
\]

Changing the adapted chart transports (GF6) by the corresponding boundary
valuation and `R`-preserving Poisson action; the invariant datum is the zero
orbit, not the displayed Laurent representative.  In the weighted family
the generic and exceptional representatives are respectively the universal
lines `Theta` and `Theta_exc` from the all-degree descent theorem.  The
quadratic or `XQ` shear kills them, proving that (GF4) vanishes throughout the
admissible family.

Thus the cotangent graph separates the dimension drop into a universal
restriction and one precise polynomial issue: vanishing of graph flux.

## 3. The natural graph centralizer is exactly the transported `A_2`

The same graph has an immediate Weyl model, but it gives a structural no-go
rather than the desired endomorphism.  Use target Weyl coordinates

\[
 (s,t,r;\partial_s,\partial_t,\partial_r)
\]

and the second-class normal pair

\[
 a=\partial_s,\qquad b=s-\partial_t,qquad [a,b]=1.   \tag{12}
\]

Its centralizer is

\[
 C=k\langle r,\partial_r,\partial_t,t-\partial_s\rangle
 \simeq A_2.                                         \tag{13}
\]

Under the inverse-Jacobian endomorphism put

\[
 a'=\delta_S,\qquad b'=S-\delta_T.                   \tag{14}
\]

The four transported centralizer generators are

\[
 R,\qquad\delta_R,\qquad\delta_T,\qquad T-\delta_S. \tag{15}
\]

They form two commuting Weyl pairs.  More strongly:

> **Graph-centralizer theorem.**  For the foundational map,
> \[
> \boxed{
> \operatorname{Cent}_{A_3}(a',b')
> =k\langle R,\delta_R,\delta_T,T-\delta_S\rangle
> =\Phi_F(C)\simeq A_2.}                             \tag{16}
> \]

Here is the exhaustion proof.  First,

\[
 \ker\delta_S\cap\ker\delta_T=k[R]                 \tag{17}
\]

inside `k[x,y,z]`.  In the function field, an element killed by the two
derivations is algebraic over `k(R)`: differentiate its minimal polynomial
over `k(S,T,R)` to remove dependence on `S,T`.  On the other hand the adapted
polynomial coordinates

\[
 Q=y+xz/3,\qquad Z=z+9Q^2
\]

give

\[
 R=2x-3x^2Q,qquad k(x,y,z)=k(R,x,Z).                \tag{18}
\]

Thus `k(R)` is relatively algebraically closed in the source function field.
Finally, the polynomial section

\[
 x=1,\quad Q=(2-R)/3,\quad Z=0                      \tag{19}
\]

shows that the polynomial intersection is exactly `k[R]`, proving (17).

Every differential operator has a unique left PBW expansion

\[
 P=\sum_\alpha c_\alpha\,
 \delta_S^{\alpha_S}\delta_T^{\alpha_T}\delta_R^{\alpha_R},
 \qquad c_\alpha\in k[x,y,z],                       \tag{20}
\]

because the inverse-Jacobian matrix is polynomially invertible.  The first
centralizer equation gives `delta_S(c_alpha)=0`.  Comparing the coefficient
of `delta^beta` in the second gives the triangular recurrence

\[
 \delta_T(c_\beta)=-(\beta_S+1)c_{\beta+e_S}.        \tag{21}
\]

At maximal differential order, the right side is zero, so (17) puts every
top coefficient in `k[R]`.  Descending in order, (21) integrates only
polynomials in `T`, with the integration constant again in `k[R]`.  Hence
all coefficients lie in `k[R,T]`, and the recurrence is exactly the normal
ordering expansion in `T-delta_S`.  This proves (16) in every differential
order.

There is an important consequence.  Any identification
`Cent(a',b') -> C` composed with `Phi_F|C` is a composition of isomorphisms,
so it is an automorphism of `A_2`.  The natural graph-normal centralizer can
never produce a non-surjective rank-two Weyl endomorphism.  Equivalently,
the polynomial Darboux coordinates on the classical graph do not extend to
new polynomial elements of this ambient Weyl centralizer.

The exact identities in (12)--(21) are checked by
[`verify_cotangent_graph_centralizer.py`](../scripts/verify_cotangent_graph_centralizer.py).

### The one-constraint normalizer gives the same no-go

Dropping the second-class partner does not enlarge the reduced algebra.  Let
`I=A_3 delta_S` be the left ideal generated by the first graph constraint and
let

\[
 N(I)=\{P\in A_3:I P\subseteq I\}
\]

be its right normalizer.  The quantum Hamiltonian reduction is `N(I)/I`
(up to the harmless opposite-algebra convention for right multiplication).

The rational inverse coordinates from the foundational map are

\[
 \tau=y+1/x,\qquad \rho=2/x,
\]

and they satisfy

\[
 T=\rho+4\tau-3R\tau^2,qquad
 2S=\tau^2+\rho\tau/2-R\tau^3.                       \tag{22}
\]

Thus

\[
 k(x,y,z)=k(T,R,\tau),qquad \ker_{k(x,y,z)}\delta_S=k(T,R).
\]

The polynomial intersection is `k[T,R]`.  Indeed, for any prime divisor in
the `(T,R)`-plane, choose `tau` generically and put
`rho=T-4tau+3R tau^2`; away from the single equation `rho=0`, the formulas
`x=2/rho`, `y=tau-1/x` and the equation for `R` give a point over its generic
point.  Hence every target prime divisor has a dominating source component,
so a denominator in `k(T,R)` cannot disappear after pullback.

In the left PBW expansion (20), reduction modulo `I` discards precisely the
terms containing a positive power of `delta_S`.  The normalizer condition is
`[delta_S,P] in I`, so the surviving coefficients lie in `ker delta_S`.
Therefore

\[
 \boxed{
 N(A_3\delta_S)/(A_3\delta_S)
 =k\langle T,R,\delta_T,\delta_R\rangle
 =\Phi_F(A_2).}                                      \tag{23}
\]

The most immediate `G_a`/moment-map reduction is consequently exhausted as
well: it returns an isomorphic transported copy, not a non-surjective
endomorphism of a fixed `A_2`.

## 4. A necessary semiconjugacy for position-type pairs

The first bounded ansatz takes

\[
 h=h(X_1,X_2,X_3),\qquad
 L=\sum_i v_i(X)D_i+c(X),\qquad [L,h]=1.
\]

The principal symbol of `L` is nonzero.  Normal ordering shows that every
nonzero expression

\[
 \sum_{j=0}^m a_j(h)L^j,\qquad a_m\ne0,
\]

has differential order `m`.  Consequently

\[
 k\langle h,L\rangle\cap k[X_1,X_2,X_3]=k[h].       \tag{1}
\]

If this Weyl subalgebra is invariant, then its multiplication generator must
satisfy

\[
 h(F)=g(h)                                                \tag{2}
\]

for a univariate polynomial `g`.  In particular,

\[
 dh\wedge d(h\circ F)=0.                                 \tag{3}
\]

Equation (3) is a cheap necessary condition independent of the unknown
conjugate operator `L`.

## 5. Exact bounded no-go result

Let `h` be the general polynomial of total degree at most four, with its
constant term removed.  Expanding (3) gives `4160` distinct homogeneous
quadratic equations in the `34` coefficients of `h`.

Over `F_32003`, a Gröbner basis has dimension zero, vector-space dimension
`35`, and size `595`.  Reduction by that basis puts the fourth power of every
coefficient in the ideal.  Hence the radical is the irrelevant maximal ideal
and the projective solution scheme is empty.

Because the equations are integral and homogeneous, a nonzero solution over
characteristic zero could be scaled to have integral coefficients at a place
over `32003`, with at least one unit coefficient.  Reduction would give a
nonzero projective solution over the algebraic closure of `F_32003`, a
contradiction.  Therefore:

\[
 \boxed{\text{No nonconstant }h\text{ of degree at most four can occur in an
 invariant position/first-order Weyl pair.}}              \tag{4}
\]

The three immediate image pairs do not evade the test.  If
`delta_i(F_j)=delta_ij`, then `k<F_i,delta_i>` is a Weyl subalgebra, but its
invariance would require `F_i o F` to lie in `k[F_i]`.  The gradient wedges
are already nonzero at small integral points:

| coordinate | point | `grad(F_i) cross grad(F_i o F)` |
|---|---:|---:|
| `F_1` | `(0,0,0)` | `(0,2,0)` |
| `F_2` | `(0,0,1)` | `(0,0,-6)` |
| `F_3` | `(0,0,0)` | `(0,-4,0)` |

Thus none of the three canonical image pairs is invariant.

Run the certificate with:

```bash
.venv/bin/python scripts/search_invariant_weyl_reduction.py --max-degree 4
```

The conclusion is deliberately bounded.  It does not exclude:

* a position generator of degree at least five;
* a pair whose two generators both have positive differential order;
* an invariant rank-two centralizer not generated by a visible spectator
  pair; or
* a Hamiltonian reduction defined only after localization or completion.

## 6. Next search

The graph-centralizer theorem removes the most natural transported pair in
all differential orders, and (23) removes its one-moment normalizer
reduction.  The next finite search should therefore use the differential-order
filtration on a **non-transported** Weyl pair rather than increase the degree
of `h` blindly.

1. Take `P,Q` of differential order at most one and bounded Bernstein degree.
2. Solve `[P,Q]=1` at principal and subprincipal order.
3. Exclude pairs polynomially generated by `(a',b')`, then impose that the
   principal symbols span an invariant symplectic two-plane under the
   cotangent symbol map of `Phi_F`.
4. For surviving components, test centralizer closure and whether the
   associated graded centralizer is a polynomial symplectic four-fold.
5. Only then solve exact invariance in the Weyl algebra and transport the
   known three-point collision to prove that non-surjectivity survives.

In parallel, the graph constraints are second-class: their commutator is a
unit.  They therefore admit no Hamiltonian-normalizer quotient.  A quantum
version of the intrinsic graph which retains the larger classical graph
coordinate ring must use a module/corner or Weyl bimodule, not the ambient
centralizer (16).  The existing formal-local construction already produces
such a rank-one formal Weyl bimodule.  The sharp global question is whether
it descends across `X=0` with polynomial coefficients.  Its pole filtration
should be compared against the same four-residue line and the existing
`hbar^5` obstruction.
