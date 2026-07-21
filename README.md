# A three-dimensional Keller map with a rational collision

This repository studies the polynomial map below and its marked-root
geometry. The stable core is C01--C04; later work is separated into
experimental areas.

## Explicit map

Let `u=1+xy` and define `F:A^3 -> A^3` by

\[
F(x,y,z)=\left(
u^3z+y^2u(4+3xy),
y+3xu^2z+3xy^2(4+3xy),
2x-3x^2y-x^3z
\right).
\]

The coefficients are rational, so the map is defined over every
characteristic-zero field.

## Determinant certificate

On `x!=0`, set

\[
t=y+1/x,\qquad r=2/x,\qquad c=F_3(x,y,z).
\]

Writing the first two target coordinates as `(a,b)` gives

\[
a=t^2+rt/2-ct^3,\qquad b=r+4t-3ct^2.
\]

The two coordinate Jacobians are

\[
\det {\partial(a,b,c)\over\partial(t,r,c)}={r\over2},
\qquad
\det {\partial(t,r,c)\over\partial(x,y,z)}=-2x.
\]

Since `r=2/x`, their product is `-2`. Both sides are polynomials, hence

\[
\boxed{\det DF=-2}
\]

on all of `A^3`.

## Explicit collision

Exact substitution gives

\[
F(0,0,-1/4)=F(1,-3/2,13/2)=F(-1,3/2,13/2)=(-1/4,0,0).
\]

The three source points are distinct. Thus `F` is an everywhere-etale,
noninjective polynomial map. Appending identity coordinates gives the same
phenomenon in every dimension at least three.

## Geometric interpretation

For a target `(a,b,c)`, form the binary cubic

\[
Q_{a,b,c}(U,V)=cU^3-2U^2V+bUV^2-2aV^3.
\]

The affine source is isomorphic to the incidence space consisting of a cubic
together with a marked simple projective root. Under this isomorphism, `F`
forgets the marked root. Consequently:

- the generic fiber has three points;
- the fiber cardinalities are `3`, `1`, and `0` for root types
  `(1,1,1)`, `(2,1)`, and `(3)`;
- the generic monodromy group is `S_3`;
- the image and nonproperness set are controlled by the cubic discriminant.

The weighted C04 construction replaces the cubic by
`H(W)-BCW+cAC^2` and retains the same normalized marked-root mechanism.

## Mathematical status

C01–C04 are the stable core. See [STATUS.md](STATUS.md) for detailed
evidence levels.

## Deeper work

- [Core paper draft](papers/core-counterexample/main.tex)
- [C01 exact proof](verified/FOUNDATIONAL_GEOMETRY.md)
- [C02 marked-root model](verified/MARKED_ROOT_MODEL.md)
- [C03 exact image and nonproperness theorem](verified/IMAGE_AND_NONPROPERNESS.md)
- [C04 weighted theorem](verified/WEIGHTED_SEED_THEOREM.md)
- [Cancellation construction](experimental/cancellation/CONSTRUCTION.md)
- [Mathematical status and external review](STATUS.md)
- [Reproduction commands](REPRODUCE.md)

## Further results and archives

- [Standalone discriminant-pencils paper](papers/discriminant-pencils/main.tex)
- [Adversarial audit of C05](papers/discriminant-pencils/AUDIT.md)
- [Experimental geometry, C05--C16](experimental/geometry/README.md)
- Cancellation programme: [arithmetic](experimental/cancellation/ARITHMETIC.md),
  [boundary geometry](experimental/cancellation/BOUNDARY_GEOMETRY.md),
  [rigidity within the current ansatz](experimental/cancellation/RIGIDITY.md),
  and [open problems](experimental/cancellation/OPEN_PROBLEMS.md)
- [Archived notes and superseded derivations](archive/README.md)
