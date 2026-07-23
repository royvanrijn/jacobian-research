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

## 2. A necessary semiconjugacy for position-type pairs

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

## 3. Exact bounded no-go result

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

## 4. Next search

The next finite search should use the differential-order filtration rather
than increase the degree of `h` blindly.

1. Take `P,Q` of differential order at most one and bounded Bernstein degree.
2. Solve `[P,Q]=1` at principal and subprincipal order.
3. Impose that the principal symbols span an invariant symplectic two-plane
   under the cotangent symbol map of `Phi_F`.
4. For surviving components, test centralizer closure and whether the
   associated graded centralizer is a polynomial symplectic four-fold.
5. Only then solve exact invariance in the Weyl algebra and transport the
   known three-point collision to prove that non-surjectivity survives.

In parallel, the classical rank-two completion gives a better constrained
quantization problem than a general invariant-pair search.  Its Ore/Darboux
localization isolates the global issue to polynomial pole descent across
`X=0`; any invariant reduction found in `A_3` should be compared against that
same four-residue line and the existing `hbar^5` obstruction.
