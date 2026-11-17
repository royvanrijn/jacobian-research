# Complete orthogonal block frames: minimum36, next value47

KCUB24, 15 November 2026. Let K be the fixed integral twelve-variable cubic
of [KCUB3](PROOF12.md), with36 nonlinear monomial occurrences. Consider
linear source matrices L with the following exact scope:

* L acts by an arbitrary orthogonal transformation of
  Q=2*x2^2-x11*x7 on the three coordinates (x2,x7,x11);
* it fixes x2+x10, adjusting x10 accordingly;
* all other eight coordinates are fixed.

Then L^-1*K(LX) has at least36 nonlinear occurrences. Equality occurs
exactly for the monomial matrices described below. There are no members
with37,...,46 occurrences. The next count47 is attained.

This completes a full algebraic subgroup, including both components and
all boundary charts, rather than a bounded list of transformations. It is
not a bound for the full affine orbit, transformations mixing the remaining
eight coordinates, arbitrary nonlinear changes or all cubic maps. No new
positive counterexample or dimension/term record is claimed. Current bounds
remain ordinary cubic12/36, cubic-linear89, raw3D13, normalized3D[10,48],
global3D degree[5,7], and exact cubic-component degree minimum9. The strict
census still has62 profiles excluded and12 open; the full goal stays open.

## Complete parameterization of the group

Put b=x2, d=x7, f=x11 and e=x10. Represent the quadratic by
S=[[f,b],[b,d/2]], whose determinant is -Q/2. For
M=[[alpha,beta],[gamma,delta]] in SL2 and epsilon in {1,-1}, set
S'=epsilon*M*S*M^T and e'=b+e-b'. This defines R_epsilon(M)=L.
Explicitly,

```
f'=epsilon*(alpha^2*f+2*alpha*beta*b+beta^2*d/2),
b'=epsilon*(alpha*gamma*f+(alpha*delta+beta*gamma)*b+beta*delta*d/2),
d'=epsilon*(2*gamma^2*f+4*gamma*delta*b+delta^2*d),
e'=b+e-b'.
```

These give the entire stated orthogonal group over C. Indeed, an orthogonal
transformation preserves the smooth projective conic detS=0. Parametrizing
rank-one symmetric matrices identifies this conic with P1. Its induced
projective automorphism is Mobius, hence comes from some M in SL2 after
scalar normalization. The two linear actions differ by a scalar; preserving
Q forces its square to be1, giving epsilon. Extension to C proves the lower
bound for every characteristic-zero coefficient field, and all minima below
are already attained over Q. The conic argument is written, not formalized.

Right multiplication of M by D=diag(lambda,lambda^-1) changes L by a right
diagonal matrix: f scales bylambda^2, d bylambda^-2, and b,e stay fixed.
Such diagonal conjugation preserves every coefficient zero and support count.
Consequently two charts suffice for each sign:

```
M_open(t,u)=[[1+t*u,t],[u,1]],
M_boundary(a)=[[a,1],[-1,0]].
```

For delta!=0, M*diag(delta,delta^-1)=M_open(beta/delta,gamma*delta).
For delta=0, beta!=0 and
M*diag(beta,beta^-1)=M_boundary(alpha*beta).
This covers all M without dropping a denominator-zero stratum. The exact
checker reconstructs all four resulting polynomial12-by12 matrices and
their inverses, both preserved forms and every coefficient of L^-1*K(LX).

## Coefficient factors and a finite universal count

The generic support counts in the plus/minus open charts are126, and in
the plus/minus boundary charts72 and60. Every coefficient is a nonzero
rational scalar times a product from the finite factor list in its model.
The checker divides each independently reconstructed coefficient by those
listed polynomials until the remainder is a nonzero rational constant.
It then reconstructs the multiplicity of every factor set. Thus the support
at every complex parameter value is exactly the sum of the corresponding
nonzero-product indicators; no numerical root sampling is used.

For the open charts write v=1+t*u and w=1/2+t*u. Since v-w=1/2,
they cannot both vanish. Counting only factor groups whose factors are
forced nonzero gives the following lower bounds. Omitting other groups
can only decrease the count.

| Open stratum | Plus sign | Minus sign |
| --- | ---: | ---: |
|t,u,v,w all nonzero|88|88|
|t,u nonzero, v=0|59|49|
|t,u nonzero, w=0|62|62|
|t=0,u nonzero|49|59|
|u=0,t nonzero|50|59|
|t=u=0 (exact count)|36|47|

The last-axis row uses the exact univariate specializations at u=0:
each selected factor becomes a nonzero constant or a nonzero monomial in t.
This accounts for all parameter values, including every coincident factor
root. The base origin is independently reconstructed as K in the plus chart.

The boundary charts have

| Boundary stratum | Plus sign | Minus sign |
| --- | ---: | ---: |
|a nonzero (lower bound)|59|50|
|a=0 (exact count)|47|36|

These rows again count only coefficients with factors forced nonzero.
The plus boundary also has38 parameter-independent occurrences, but the
sharper split above is needed for the gap statement.

Thus36 is the minimum in the whole group. All other chart points have
at least47 occurrences, and47 occurs in the minus open origin and the
plus boundary origin. The only36 representatives are the plus open origin
and the minus boundary origin. Restoring the right diagonal factor, the
first is a diagonal matrix; the second sends

```
b -> b, f -> -d/2, d -> -2*f, e -> e,
```

followed by such a diagonal. Both are monomial source matrices, so their
conjugates preserve the original36-term support. This proves the exact
minimum and sharp next count. It does not characterize minima outside this
defined block group or prove any global37-to46 gap for cubic Keller maps.

## Independent replay and retained failed route

Run

```
.venv/bin/python research/jacobian-followup/cubic-reduction/verify_orthogonal_surface12.py
```

The checker uses only standard-library Fraction arithmetic and a guarded
AST parser, reusing the pinned sparse arithmetic in verify.py. It constructs
the matrices from the displayed formulas, composes the pinned K, checks
both inverse identities and all coefficient dictionaries, verifies every
factor product and the complete indicator counts. SymPy factorization in
the primary producer is not trusted as a proof. The written conic/chart
arguments are separate from these exact arithmetic checks. The determinant
and collision of every conjugate follow by chain rule from pinned KCUB3;
this bound does not need a new determinant expansion in all parameters.

- [Independent receipt](orthogonal-surface12/validation.json)
- [Open plus model](orthogonal-surface12/open/model.json)
- [Boundary plus model](orthogonal-surface12/boundary/model.json)
- [Open minus model](orthogonal-surface12/open-minus/model.json)
- [Boundary minus model](orthogonal-surface12/boundary-minus/model.json)
- [Primary producer](orthogonal_surface12.py)

The first calculation of this turn, `base-quadratic-conjugacy/`, was
redundant: KCUB19 already excludes arbitrary-degree one-axis nonlinear
source shears, including all29 combined base-axis atoms tested there.
Its primary bases and an initial Singular fraction-renderer error are
retained, not promoted as a new theorem. The renderer was fixed by writing
rational coefficients before integer monomial powers, with the original
input/log/source snapshot preserved. Do not rerun these one-axis models.
The new mathematical progress is the complete signed orthogonal block
cover above, beyond the previously saved one-parameter line.

A further support reduction must leave this group: e.g. change the
quadratic form or the fixed sum, mix additional weight blocks, use genuinely
coordinated nonlinear source/target operations, or change the construction.
No full-affine or global optimum follows from the present obstruction.
