# Complete weighted cubic carrier realization (KCUB20)

Fix the integral12-variable cubic K of [KCUB3](PROOF12.md), and its nine
carrier outputs C=(K4,...,K12). In the original triangular order4,...,12,
allow **every** weight-compatible monomial of degrees2 and3 in preceding
coordinates in each carrier output, keeping linear partI. There are135
independent coefficients. For the resulting carrier map C_z, define the
triangular source automorphism E_z by C(E_z(X))=C_z(X), fixing x1,x2,x3,
and put G_z=K E_z.

The complete cubic family is the eight-parameter family below. Its minimum
number of displayed nonlinear occurrences is36, attained only by diagonal
conjugates of K. It has no37- or38-term member. Every member has twelve
linearly independent nonlinear output parts, excluding an affine-linear
output component in any affine frame. The cubic-output rank is at least5;
the rank5 locus is linearly conjugate to K. Thus the standard homogeneous
companion recipe has minimum dimension18 within this family and supplies
no new tensor on that minimum locus.

These are scoped construction limits. They do not prove a global cubic
term or dimension minimum, exclude nonlinear graph destabilizations, or
bound all carrier orders, changes of base variables or other constructions.
No smaller counterexample was found. [KCUB19](CARRIER_REALIZATION_LIMIT.md)
remains the earlier47-atom subfamily theorem and the separate all-degree
one-axis nonlinear source-shear obstruction.

## Complete model and trusted base

The weights are

```
(1,-1,-2,0,-3,2,-1,-2,-2,-1,-1,0).
```

For i=4,...,12 enumerate all exponent vectors e with |e|=2 or3,
e_j=0 for j>=i, and sum_j weight_j e_j=weight_i. The independent checker
re-enumerates exactly135 atoms. The original12 atoms retain their old
coefficient as a multiplier of z_j; the others have coefficient z_j.
The complete labeled table is in `carrier-reparam-full/input.json`.
The first47 labels agree with KCUB19; labels48,...,135 are the88 previously
missing carrier-dependent cubic atoms.

The old and new carrier changes are triangular over the three base
variables with unit diagonal. Hence every E_z is a polynomial automorphism
with E_z(0)=0 and JE_z(0)=I. Every G_z inherits det DG_z=1 and the collision
through E_z^-1 from KCUB3. Those base facts are imported from its pinned
witness and validation, not re-proved here. Cubic degree requires exactly
the vanishing of540 distinct coefficient equations of source degree>3.
Before solving them there are184 nonlinear slots of degree2 or3.

## Exact complete parameterization

The certificate expresses127 necessary relations in the ideal of the540
high-degree equations. They consist of102 individual z_j=0 for j outside

```
A={1,2,3,4,5,6,7,8,9,10,11,12,13,17,18,19,22,26,27,28,30,
   34,35,38,39,42,45,79,80,83,96,117,118}
```

and the following25 relations:

```
z3-1, z1-z2, z7-z9, z8-z9, z10-z12,
z2*z9-1, z5*z12-z1, z4-z9*z12, z11-z4-3*z22,
z5*z42+z22, z2*z17-2*z13, z2*z18-z5*z22,
z5*z38-3*z5*z96-z17,
z12*z30-z96, z12*z34-z117, z12*z35-z118,
z2*z80-z27, z9*z28-z118, z9*z26-z79,
z12*z79-z12*z117+3*z28, z9*z19+3*z34,
z6*z9-z9-z19, z5*z39+z96, z5*z45-z80, z83-z118.
```

All127 identities are multiplied independently, using410 multiplier
monomials, maximum multiplier degree5. No generic parameter denominator
is discarded. The product identities force z9,z2,z12,z5 nonzero. Set

```
p=z9, q=z12, s=z17, t=z22, u=z117, v=z118, w=z80, z=z96.
```

Then pq!=0 and the remaining coefficients are forced to be

```
z1=z2=1/p, z3=1, z4=p*q, z5=1/(p*q),
z6=1-3*u/(p^2*q), z7=z8=z9=p, z10=z12=q,
z11=p*q+3*t, z13=s/(2*p), z17=s, z18=t/q,
z19=-3*u/(p*q), z22=t, z26=(u-3*v/(p*q))/p,
z27=w/p, z28=v/p, z30=z/q, z34=u/q, z35=v/q,
z38=p*q*s+3*z, z39=-p*q*z, z42=-p*q*t, z45=p*q*w,
z79=u-3*v/(p*q), z80=w, z83=v, z96=z, z117=u, z118=v;
all other z_j=0.
```

Conversely, exact substitution makes every high-degree coefficient zero.
This proves completeness, including every zero/nonzero free-coefficient
stratum. The elementary elimination divides only by p,q whose inverses
are forced by the product identities. Equivalently, the coefficient ideal
has quotient Q[p,p^-1,q,q^-1,s,t,u,v,w,z]: the displayed necessary relations
define that Laurent polynomial ring, and every original high generator
vanishes there. This statement follows from the explicit identities and
elimination, not the primary basis dimension alone.

## Normalize scaling and count every occurrence

Let

```
D=diag(1,1,1,p,1,1/(p*q),p*q,1,1/p,1/p,1/(p*q),1/(p*q)),
S=s/p, T=t/(p*q), U=u/(p^2*q), V=v/(p^3*q^2),
W=w/p^2, Z=z/(p^2*q).
```

Coefficientwise verification gives

```
G(p,q,s,t,u,v,w,z)(X)=D^-1 G(1,1,S,T,U,V,W,Z)(D X).
```

Diagonal changes preserve support, so it suffices to use p=q=1. Below the
six normalized parameters are again denoted s,t,u,v,w,z. The generic map
has64 nonlinear occurrences. Up to nonzero rational scalars, their
coefficient inventory is

| Coefficient factor | Multiplicity |
| --- | ---: |
|1|31|
|u-1|1|
|t-u/3+v+1/3|1|
|t+w|1|
|t+v|1|
|u|4|
|v-w|1|
|s+3/2|1|
|s+3z|2|
|z|5|
|s|2|
|u-1/3|1|
|u-3v|2|
|w|3|
|v|4|
|t|3|
|t+1/3|1|

Every factor and multiplicity is independently reconstructed. Divide the
nonconstant factors into A(s,z) and B(t,u,v,w). A has weights
[s+3/2]+2[s+3z]+5[z]+2[s], with brackets meaning nonzero indicators.
It is1 at s=z=0. If z!=0 it is at least5; if z=0,s!=0 it is at least4.

B is4 at t=u=v=w=0 and is at least7 otherwise:

- If u!=0 and v!=0, weights4[u]+4[v] give8. If u!=0,v=0,
  weights4[u]+2[u-3v] give6 and at least one of u-1,u-1/3 is nonzero.
- If u=0,v!=0, weights4[v]+2[u-3v] give6 and u-1,u-1/3 add2.
- If u=v=0,w!=0, weights3[w]+[v-w] and u-1,u-1/3 give6;
  t+w and t+v=t cannot both vanish, giving at least1 more.
- If u=v=w=0,t!=0, weights3[t]+[t+w]+[t+v] and
  u-1,u-1/3 give7.

Therefore total support31+A+B is at least36, with equality exactly when
all six normalized parameters vanish. Any other point has at least39
occurrences. The39 bound is attained, for example, at s=-3/2 and the other
five zero. At equality the map is exactly the diagonal conjugate of K.
This is a count in the displayed carrier realization; arbitrary further
coordinate changes are outside the support theorem.

The primary support<=35 affine-arithmetic query is retained, with UNSAT
in its exact p=q=1 model. It is not a proof premise; the factor count and
inequalities above provide the independent argument.

## No affine-linear output component

Consider the twelve nonlinear parts G_i-x_i as polynomials. We prove their
linear independence over the coefficient field for every family member.
The weights separate the following coefficient-minor blocks:

| Output rows | Monomial columns | Determinant |
| --- | --- | --- |
|1|x1^2*x2|-3|
|5|x2^3|1|
|6|x1^2|p*q|
|4,12|x1^2*x3; x1*x10|q/p|
|2,7,10,11|x4*x10; x1*x3; x1*x8; x1*x9|12|

All are nonzero. For the remaining rows3,8,9, write A,B,C for the
coefficients of x2^2,x2*x7,x7^2 in G8. They are

```
A=1-3u/(p^2q), B=u/p-3v/(p^2q), C=v/p,
A+3B/(pq)+9C/(p^2q^2)=1.
```

At least one is nonzero. The3-by3 minor using columns
x1*x4*x5, x1*x5, and respectively x2^2,x2*x7,x7^2 equals
-16p^2 times A,B,C. Thus this weight block also has full rank.
All twelve nonlinear parts are independent.

Consequently no nonzero constant linear combination of outputs has degree
at most1. Invertible affine input changes preserve the degree of any
polynomial, and invertible affine output changes take exactly such linear
combinations plus constants. Thus no affine frame produces an affine-linear
component. This excludes deletion of a linear component, not nonlinear
carrier graph restrictions, general polynomial transformations or global
eleven-dimensional counterexamples.

## Cubic rank and the homogeneous companion

In normalized coordinates, five cubic monomial pivots in rows1,...,5 have
determinant-288, so the homogeneous cubic output part has rank at least5.
Additional same-weight cubic minors are8s, u, v, 12w, 12z. They force
rank>5 if any of s,u,v,w,z is nonzero. Conversely, when all five vanish,
only rows1,...,5 and12 have cubic terms, and the last is t times row4.
Thus the exact rank5 locus is

```
s=u=v=w=z=0, pq!=0, t arbitrary.
```

On p=q=1 in this locus put

```
S(t)=I+t E_(10,7)+t E_(12,4)       (one-based indices).
```

The entire polynomial map, including quadratic parts, satisfies
G_t=S(t) K S(t)^-1. This identity is checked independently. Together with
the diagonal normalization, every rank5 member is linearly conjugate to K.

The [standard companion construction](POWER89.md) for n variables and
cubic-output rank r uses n+r+1 variables. For this family its minimum is
therefore12+5+1=18; all points attaining it supply a linearly equivalent
companion tensor to the existing one after transporting the cubic-output
basis. This is not a lower bound on all homogeneous reduction methods,
all cubic-linear counterexamples or all shared-cube decompositions.
In particular it does not prove89 globally optimal. A genuinely different
companion from this family must have r>=6 and at least19 variables under
this recipe; it could still merit a new simultaneous cube decomposition.
The subsequent [KCUB21 study](CHANGED_COMPANION_PENCILS.md) tests two such
axes with exact derivative rank37 and98-cube point identities, without
improving the existing89-variable construction or closing other tensors.

## Reproduction and primary changed-order work

```
.venv/bin/python research/jacobian-followup/cubic-reduction/verify_carrier_full.py
```

[Validation](carrier-reparam-full/validation.json) uses standard-library
Fraction sparse arithmetic, guarded polynomial substitution and explicit
Laurent coefficient units. It reconstructs all135 atoms, the source
change and full composition,540 high equations and184 low slots,
127 ideal memberships (410 multiplier monomials,maxdegree5), the complete
family and its normalization, all coefficient factors, nonlinear-output
rank minors and their Bezout identity, cubic-rank minors and the rank5
full-map linear conjugacy. Written elementary elimination/counting and
companion covariance are mathematical proofs, not formalized in a proof
assistant. SymPy/Singular/Z3 are primary experiment tools only.

Ten changed triangular orders were also prepared with complete weighted
quadratic/cubic monomial supports. All primary identity-carrier and
same-weight carrier-combination queries give unit ideals; two orders
already have unit cubic-completion ideals. These finite primary results
are retained in `carrier-order-search/`. They are not independently replayed
certificates and do not exhaust all orders or all linear combinations of
all twelve outputs. Last-pivot and singleton combination charts are exactly
the separate identity-carrier cases; the rank script covers the four other
carrier charts.

A joint linear conjugacy line preserving2*x2^2-x11*x7 and x2+x10 was
screened at every rational coefficient-cancellation parameter. Its best is36
at zero; special nonzero parameters have54 or59 terms. This finite line
is retained in `orthogonal-line12/` and proves no arbitrary-frame optimum.

Current construction bounds stay ordinary cubic12 variables/36 terms,
cubic-linear89 variables, raw3D13 terms and normalized3D[10,48]. Next
change the base construction/frame/order or use a genuinely new companion
with a full determinant/collision replay. The135 original-order cubic
carrier coefficients and the rank5 companion locus are now closed routes.
