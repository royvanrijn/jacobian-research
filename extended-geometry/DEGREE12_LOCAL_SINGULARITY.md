# The generic degree-twelve local singularity

Work over a field of characteristic zero.  Let `bar(xi)` be a geometric
generic point of the exact collision stratum `E_(6,6)` in the component
`C_(3,2)`.  The two points of the normalization above it correspond to the
allocations

\[
((3,0),(0,2)),\qquad ((0,2),(3,0)).
\]

## Theorem

Etale-locally along a dense open subset of `E_(6,6)`, choose regular
parameters on the two normalization sheets so that

\[
B_+=k[[t,x,y,z]],\qquad B_-=k[[t,x',y',z']].
\]

Then the completed component ring is

\[
\boxed{
A=B_+\mathop{\times}_{D}B_-,
\qquad
D=k[[t,\epsilon,\eta]]/(\epsilon^2,\eta^2),}
\]

where the quotient maps have kernels

\[
I_+=(z,x^2,y^2),\qquad I_-=(z',x'^2,y'^2).
\]

The conductor in the completed normalization is

\[
\boxed{\mathfrak c=I_+\oplus I_-\subset B_+\oplus B_-.}
\]

After a transverse slice to `E_(6,6)`, the intersection algebra has length

\[
\boxed{\operatorname{length}k[\epsilon,\eta]/
       (\epsilon^2,\eta^2)=4.}
\]

Thus the two branches are quadratically tangent, with generic transverse
intersection multiplicity four.  They are not nodal.

## 1. The correct coefficient relation

For `M=Q^2R^3`, the normalized seed is

\[
H={-M+M'(0)W+M(0)\over D},
\qquad D=M'(1)-M'(0).
\]

The leading coefficient of `H` recovers `D`.  Equality of two normalized
seeds is therefore equality of the coefficients of their monic polynomials
`M_+,M_-` in degrees `2,...,12`.  Equivalently,

\[
M_+-M_-=\lambda W+\mu.
\]

It would be incorrect to replace this relation silently by `M_+=M_-`.

## 2. A closed length-four subsystem

First impose the stronger equation `M_+=M_-`.  Hensel separation at the two
distinct sixfold roots splits it into two independent factorization blocks

\[
Q^2=T^3,
\]

with `Q` monic cubic and `T` monic quadratic.  Translate the root and depress
the cubic:

\[
Q=Z^3+uZ+v,qquad T=Z^2+aZ+b.
\]

Coefficient comparison gives

\[
\begin{aligned}
[Z^5]&=-3a,\\
[Z^4]&=-3a^2-3b+2u,\\
[Z^3]&=-a^3-6ab+2v,\\
[Z^2]&=-3a^2b-3b^2+u^2,\\
[Z]&=-3ab^2+2uv,\\
[1]&=-b^3+v^2.
\end{aligned}
\]

Their exact Groebner basis in characteristic zero is

\[
a,\qquad v,\qquad 2u-3b,\qquad b^2.
\]

One block is therefore a root-position coordinate times one dual number.
The two roots give independent parameters `epsilon,eta` with
`epsilon^2=eta^2=0`.  The incidence equation is smooth in the two root
positions on a dense open set: at the exact witness with roots `-6/5,1`, its
two root derivatives are

\[
-{7776\over625},\qquad -{279936\over15625}.
\]

It removes one smooth root-position coordinate and no nilpotent direction.
Thus the stronger-equality subsystem has transverse algebra

\[
D_0=k[\epsilon,\eta]/(\epsilon^2,\eta^2)
\]

and length four.

## 3. The affine-difference equations add nothing

Now return to the correct relation `M_+-M_-` affine.  The two normalization
differentials have rank four.  Their image tangent spaces intersect in
dimension three, while the collision stratum has dimension one.  Hence a
transverse slice of the full off-diagonal double scheme has embedding
dimension two.

Modulo the combined image tangent space, the paired second derivatives have
rank two.  In exact common-tangent coordinates their two scalar forms can be
taken as

\[
\begin{aligned}
q_1={}&186624x_0^2+648000x_0x_1-388800x_0x_2
 +953125x_1^2-1612500x_1x_2+765000x_2^2,\\
q_2={}&(5x_1-6x_2)^2.
\end{aligned}
\]

On `5x_1-6x_2=0`,

\[
q_1=9(144x_0+125x_1)^2.
\]

After quotienting by the collision direction, these are squares of two
independent linear forms.  The associated graded transverse algebra is
therefore a quotient of

\[
k[\epsilon,\eta]/(\epsilon^2,\eta^2),
\]

and has length at most four.  But the stronger-equality algebra `D_0` is a
closed subscheme of the full double scheme and has length four.  The induced
surjection of local Artin rings has equal finite lengths, so it is an
isomorphism.  This proves

\[
D=D_0
\]

and excludes additional nilpotents and embedded components on the generic
open set.

## 4. Fiber product and conductor

The component is reduced and of finite type over a characteristic-zero
field; excellence makes its completion reduced.  At `bar(xi)` it has exactly
the two minimal primes supplied by the normalization fiber.  For any reduced
ring with exactly two minimal primes `p_+,p_-`, the standard exact sequence

\[
0\longrightarrow A\longrightarrow A/p_+\oplus A/p_-
 \longrightarrow A/(p_++p_-)\longrightarrow0
\]

identifies it with the fiber product of its two branches over their
scheme-theoretic intersection.  Sections 2--3 identify that intersection
with `D`, proving the displayed presentation of `A`.

For a fiber product of two normal rings over a common quotient, the largest
ideal of `B_+\oplus B_-` contained in the fiber product is the direct sum of
the two quotient kernels.  Hence

\[
\mathfrak c=I_+\oplus I_-.
\]

The basis `1,epsilon,eta,epsilon*eta` gives the transverse length four.

## Executable certificate

Run

```bash
python scripts/verify_degree12_branch_intersection.py
```

The script checks the two differential ranks, their tangent intersection,
the exact second quadratic forms, the sixfold-block Groebner basis, smoothness
of the incidence equation at an admissible witness, and the matching
upper/lower length bounds.

### Independent algebra audit

The decisive sixfold block also has an elementary proof independent of the
Groebner calculation.  For

\[
Q=z^3+uz+v,\qquad T=z^2+az+b,
\]

coefficient comparison in `Q^2-T^3` gives the exact coefficient ideal

\[
(c_5,c_4,c_3,c_2,c_1,c_0)=(a,v,2u-3b,b^2).
\]

Thus each separated sixfold root contributes one dual-number block, and the
two blocks give
`k[epsilon,eta]/(epsilon^2,eta^2)` of length four.  A dependency-free checker
also recovers the two common-tangent quadratic restrictions, giving the
matching upper bound.

Run

```bash
python3 scripts/audit_degree_twelve_independent.py
```

The former standalone audit narrative is retained in
[archive/geometry-support](../archive/geometry-support/README.md).
