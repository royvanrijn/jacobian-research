# Foundational Keller map and rational collision

Put `u=1+xy` and define

\[
F(x,y,z)=\left(
 u^3z+y^2u(4+3xy),
 y+3xu^2z+3xy^2(4+3xy),
 2x-3x^2y-x^3z
\right).
\]

## Theorem

Over every characteristic-zero field,

\[
\det DF=-2
\]

and

\[
F(0,0,-1/4)=F(1,-3/2,13/2)=F(-1,3/2,13/2)=(-1/4,0,0).
\]

Thus `F` is an everywhere-etale, noninjective polynomial map.  Appending
identity coordinates gives the same phenomenon in every dimension at least
three.

## Exact proof

On the dense open `x!=0`, set

\[
t=y+1/x,\qquad r=2/x,\qquad c=F_3(x,y,z),
\]

and write the first two target coordinates as `(a,b)`.  Direct substitution
gives

\[
a=t^2+\frac{rt}{2}-ct^3,
\qquad
b=r+4t-3ct^2.
\]

The two coordinate Jacobians are

\[
\det\frac{\partial(a,b,c)}{\partial(t,r,c)}=\frac r2,
\qquad
\det\frac{\partial(t,r,c)}{\partial(x,y,z)}=-2x.
\]

Since `r=2/x`, their product is `-2`.  Both sides are polynomials, so equality
on the dense open proves the determinant identity everywhere.  The collision
is exact substitution, and the three displayed source points are distinct.

This proof is intentionally confined to the foundational map.  The
nonduplicated continuation of the core chain is:

- [Normalized factorization model](NORMALIZED_FACTORIZATION_MODEL.md): the
  canonical algebraic construction and unequal-degree extension;
- [Foundational incidence construction](FOUNDATIONAL_INCIDENCE_CONSTRUCTION.md):
  projective hyperplane geometry and the exceptional `(2,1)` orbit;
- [Cubic marked-root model](MARKED_ROOT_MODEL.md): the marked-root
  interpretation and its two reconstruction charts;
- [Exact image and nonproperness](IMAGE_AND_NONPROPERNESS.md);
- [Weighted marked-root theorem](WEIGHTED_SEED_THEOREM.md).

## Reproduction

Run

```bash
make verify-minimal
make verify-core
```

The minimal target uses a dependency-free exact implementation.  The full
core target adds the normalized factorization, marked-root, and image
regressions.  External formal verification is supplied both by the separately
authored pinned [Lean certificate](LEAN_FOUNDATIONAL_MAP.md) and by the
refereed [Archive of Formal Proofs Isabelle/HOL
entry](https://isa-afp.org/entries/Jacobian_Counterexample.html).  The AFP
entry independently verifies the determinant, exact collision,
determinant-one normalization, identity-padded stabilization in every finite
dimension at least three, and nonexistence of a polynomial inverse.

The former all-in-one derivation and presentation-invariance audit are
retained in [archive/core-support](../archive/core-support/README.md).
