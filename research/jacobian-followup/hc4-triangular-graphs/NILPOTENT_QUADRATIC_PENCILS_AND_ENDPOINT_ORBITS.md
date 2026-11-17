# Nilpotent quadratic pencils have exactly the quadratic endpoint orbit

15 November 2026. Claim **HC4MYJD0**. Together with
[JC0/1](SEMISIMPLE_QUADRATIC_PENCILS_AND_ALL_SOURCE_ORBITS.md),
this separates the commonly centered quadratic-composite graph
frontier by its actual pencil: nonzero semisimple pencils forbid
unit Hessian; nilpotent pencils give exactly the old nondegenerate
quadratic source orbit. A pencil with a nonzero eigenvalue and a
nontrivial Jordan block is not classified here.

## 1. All outer degrees reduce by polynomial shears

Let q_A=x^T A x be nondegenerate in four complex variables,
q_B=x^T B x, and assume T=A^(-1)B is nilpotent. For every polynomial
H, there is a polynomial automorphism Phi with

    q_A+H(q_B)=q_A composed with Phi+H(0).                (1)

The native quadratic can first be put into a constant linear normal
form; conjugating the written shears returns (1) in the original
coordinates. There is no bound on deg H.

Here is the complete reduction, including the linear-algebra input.
The endomorphism T is self-adjoint for A. Let l be its largest
Jordan-chain length. The symmetric bilinear form
(v,w) -> A(T^(l-1)v,w) is nonzero, so characteristic zero gives
a v with A(T^(l-1)v,v)!=0. Its cyclic space

    span(v,Tv,...,T^(l-1)v)

is nondegenerate: its Gram matrix is anti-triangular with nonzero
anti-diagonal. It is T-invariant, and its A-orthogonal complement
is T-invariant as well.

One can normalize this chain so its Gram matrix is exactly
anti-diagonal with entries one. Replace v by u(T)v for a polynomial
unit u modulo t^l. First scale u(0) to normalize
A(v,T^(l-1)v). Descending through the lower pairings, the new
coefficient at each stage occurs with twice that nonzero top
pairing, and so is uniquely solvable. Higher powers do not affect
pairings already fixed. This makes A(v,T^j v)=0 for j<l-1.
Thus the block has

    q_A=sum_{i+j=l+1} x_i x_j,
    q_B=sum_{i+j=l}   x_i x_j,                          (2)

with 1<=i,j<=l; reversing the cyclic basis if needed gives this
orientation. Splitting the orthogonal complement proves an
orthogonal direct sum of these blocks. No unproved classification
of polynomial automorphisms or embeddings is used.

The partitions of four now give the following complete list,
up to constant linear equivalence:

| Jordan lengths | q_A | q_B |
|---|---|---|
| 1+1+1+1 | x^2+y^2+z^2+w^2 | 0 |
| 2+1+1 | 2xy+z^2+w^2 | x^2 |
| 2+2 | 2xy+2zw | x^2+z^2 |
| 3+1 | 2xy+z^2+w^2 | 2xz |
| 4 | 2xw+2yz | 2xz+y^2 |

Write H0=H(0). The zero block needs only the identity.
For length 2+1+1 put

    R=(H(x^2)-H0)/(2x),        Phi=(x,y+R,z,w).

For length 2+2 put

    R=(H(x^2+z^2)-H(z^2))/(2x),
    S=(H(z^2)-H0)/(2z),
    Phi=(x,y+R,z,w+S).

For length 3+1 put

    R=(H(2xz)-H0)/(2x),        Phi=(x,y+R,z,w).

In every case the displayed quotient is polynomial: its numerator
vanishes at the corresponding x=0 or z=0. The shifted variables
do not occur in the shifts, so the inverse subtracts those shifts.
Substitution in q_A gives (1) directly.

For the length-four block put

    R=(H(2xz+y^2)-H(y^2))/(2x),
    S=(H(y^2)-H0)/(2y),
    Phi=(x,y,z+S,w+R).                                  (3)

Again R and S are polynomials. This is triangular, even though
R can depend on the old z. Its inverse is

    x=X, y=Y, z=Z-S(Y),
    w=W-R(X,Y,Z-S(Y)).                                  (4)

Now

    q_A composed with Phi
      =2xw+2yz+2xR+2yS
      =q_A+H(2xz+y^2)-H0.

This proves (1) in the last block and hence in every case. The
normalizing source map is a composition of constant linear maps
and explicit polynomial shears; no formal or rational inverse
is substituted for a polynomial automorphism.

## 2. Exact full-source-orbit and graph-image consequence

Equation (1) identifies the full polynomial source-automorphism
orbit of q_A+H(q_B) with the full orbit of q_A+H(0). It does not
make Hessians invariant under the normalizing shear. Instead the
sets of actual source potentials are equal because composing
with a fixed polynomial automorphism permutes all source
automorphisms.

By the retained
[HG2 quadratic-orbit theorem](GLOBAL_VALUE_ORBITS.md#3-the-two-known-endpoints-are-full-source-orbits-hc4myhg2),
the universal gradient-inverse question for unit-Hessian members
of this orbit is exactly the singleton-zero plane Keller endpoint.
The point is one reduced full critical point, preserved by every
source automorphism. Its existence does not prove plane Keller
properness or invertibility.

Apply this to any fixed actual graph

    S:X=P0(q_B(u)) in A5,
    Q(u,X)=q_A(u)+alpha X^2+beta X+gamma,
    alpha!=0, det A!=0,                                 (5)

with nilpotent A^(-1)B and arbitrary P0. Set
H=alpha P0^2+beta P0+gamma. Every closed parameterization of S
factors through its graph map by a polynomial source automorphism.
Conversely every such automorphism gives a closed parameterization.
Therefore this one graph image already has exactly the old
nondegenerate quadratic source orbit, up to its harmless additive
critical value H(0). This is an equality of orbits and a conditional
endpoint realization, not an asserted solution of that endpoint.

The equality transfers through every certified polynomial
quadratic similitude by retained
[HP0](QUADRATIC_ISOMETRY_SATURATION.md#1-saturation-under-every-polynomial-quadratic-similitude-hc4myhp0).
A nonzero scalar multiple of a complex nondegenerate quadratic
has the same source orbit after a constant scaling, up to an
additive constant. The
[HM0 literal normalization](AMBIENT_QUADRIC_GRAPHS.md#1-a-literal-ambient-automorphism-hc4myhm0)
then transfers these actual images to the original Meng--Yang
potential. Neither step assumes nonlinear gradient conjugacy.

## 3. Exact boundaries and the remaining spectral gate

The n=4 restriction belongs to the explicit complete block list;
no all-dimensional triangular normal-form theorem is asserted.
The forms must have a common center. The only eigenvalue of T
must be zero. Merely being nonsemisimple is not enough for (1).

Together JC0/1 and JD0 give, for nonzero B, nonlinear H and
commonly centered four-variable quadratic-composite graphs:

- Diagonalizable T: no source parameterization has unit Hessian.
- Nilpotent T: the full source orbit is exactly the existing
  singleton-zero endpoint.
- Nondiagonalizable T with a nonzero eigenvalue: unresolved here.

Constant graph P0 already gives the same hyperplane/quadratic
endpoint; independent shifts and arbitrary images are separate
gates. The last bullet is a structural pencil question, not a
request to scan larger graph or parameterization degrees.

The length-four native potential illustrates why (1) must not
be presented as Hessian invariance. For H(t)=t^2 its native
Hessian determinant is

    16*(1+8xy-8x^2*(2xz+y^2)),

which is nonconstant, although composing with the explicit inverse
of (3) gives the quadratic potential with constant Hessian.
The checker verifies this sign and coefficient along with the
polynomial division, inverse and normal-form identities for all
four nonzero block types. The written factor-divisibility and
chain argument, not those finite outer degrees, prove the
unbounded statement.

No new unproved external statement is imported. The proof uses
elementary self-adjoint linear algebra, exact polynomial shears,
and the retained HG2, HP0 and HM0 results. The unrestricted HC4
and plane JC2 problems remain open.
