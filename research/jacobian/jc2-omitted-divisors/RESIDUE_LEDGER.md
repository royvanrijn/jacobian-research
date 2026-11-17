# Residue degrees and every residual r7 fibre-length pattern

14 November 2026. PF2OMIT6 refines the three residual immersed-r7 rows
of [PF2OMIT5](R7_INFINITY_COMPLETION.md). The other ramified residue
degree is now proved to be one. Its target equation remains unknown.
The earlier residual nulls are retained checkpoints, superseded by this
[current ledger](residue_ledger.json), not assertions of absence.

## 1. Residue ramification is part of the exceptional ambient length

Use the finite plane normalization and global boundary hypotheses of
[PF2OMIT2](GLOBAL_LEDGER.md). Each omitted prime D is unibranch, with
normalization A1. Its image C_D has normalization A1. The induced finite
polynomial residue map g_D:A1->A1 has degree f_D. Let q_D(a) be its
local degree at the normalization point a of D; unibranchness identifies
a with a unique point of D. The generic transverse normal index is e_D,
and the ambient local fibre length at a is ell_a.

Conservation of local degree gives

    ell_a >= e_D*q_D(a).                              (R1)

Indeed choose small neighborhoods isolating a in the finite fibre.
After restricting to a nearby smooth target-image value, the residue
point splits into q_D(a) points of D near a. Each has ambient length
e_D. The local finite flat degree is conserved and any other nearby
sheets add a nonnegative amount. This works for a singular normal N
as well: its finite flat algebra has constant local cluster rank.
Only the reduced support is used to isolate the q_D(a) residue points.

Define eta_a=ell_a-e_D*q_D(a)>=0. Riemann--Hurwitz for the polynomial
residue map, including its total ramification at infinity, gives

    sum_(finite a)(q_D(a)-1)=f_D-1,
    E_D=e_D*(f_D-1)+sum_a eta_a.                       (R2)

Consequently the prime-indexed Euler formula becomes the exact identity

    d-1=sum_D e_D*f_D + sum_(D,a) eta_a.               (R3)

This reconciles the ambient finite-normalization ledger with the
already retained [residue-normalized Orevkov argument](../../foundations/plane-jc/OREVKOV_RESIDUE_DEGREE_BUDGET.md).
It is not a new attribution of that identity. In particular f is not
inserted into the old formula while retaining its old E_D unchanged.
Equations (R1)--(R2) explain exactly how E_D changes.

## 2. Zero excess also fixes image immersion

We import Orevkov, [On three-sheeted polynomial mappings of C2](https://www.math.univ-toulouse.fr/~orevkov/jc86.pdf),
Lemma5.2, with its general quantifiers, not merely its degree-three
application. At a finite boundary point where local multiplicity equals
the generic multiplicity of its boundary prime, the induced image germ
of the boundary normalization is a nonsingular embedding.

The model comparison is essential. In Orevkov Section4, the curves
contracted over finite target points are exactly the connected proper
fibres of the resolved map; the target is still the original affine A2.
Their contraction is the affine analytic Stein factor. Its normalization
is our N, by uniqueness of finite normalization. Thus his finite local
multiplicity is our ambient flat fibre length ell_a, and his generic
multiplicity is e_D. The one-point compactification at infinity plays
no role in this identification at finite points. His proof uses the
linear resolution chain and the complement of its end curve, not a
source-finality assertion. The primary pages7--9 were inspected and
[pinned](sources/residue_manifest.json); the lemma is imported, not
independently reproved by our checker.

It follows that E_D=0 implies f_D=1 by (R2), and every finite image
branch is immersive by Lemma5.2. Conversely f_D=1 and an everywhere
immersive image normalization give the local cyclic cover of PF2OMIT2,
so ell_a=e_D and E_D=0. More locally, under f_D=1 an immersive branch
has no excess, while every nonimmersive parameter has positive excess.
Different omitted primes and different branches of the same prime
remain disjoint at their normalization points.

This is specific to a globally realized Keller boundary. The abstract
double cover h^2=u^3-v^2 has ambient length2 at its cusp, equal to its
generic index2. It is not a counterexample: that boundary germ does not
occur as this distinguished boundary of a plane Keller normalization.

## 3. All residue degrees in the residual r7 table are one

Assume a degree-six plane Keller map contains an immersed e2,f1 prime R
with birational degree(6,10), gap21 image C. Additional ramification is
allowed. PF2OMIT5 excludes the entire-single-prime alternative and gives
exactly one other ramified prime S, with image B different from C.
There is at most one unramified omitted prime D, with image H if present.
Put E_S for its ambient excess and t for the number of these unramified
primes. The complete possibilities are

| Row | e_S | E_S | t | f_S forced by (R2) | Image normalization of B |
|---|---:|---:|---:|---:|---|
| A | 2 | 0 | 1 | 1 | immersive everywhere |
| B | 2 | 1 | 0 | 1 | exactly one nonimmersive parameter |
| C | 3 | 0 | 0 | 1 | immersive everywhere |

For the middle row, 2*(f_S-1)<=1 forces f_S=1. In every row E_R=0
and f_R=1; in row A the omitted unramified D has e_D=f_D=1 and an
immersive image. Thus every omitted residue degree is now known.

In row B there is exactly one point a* on S with ell_(a*)=3; every
other point has length2. Lemma5.2 makes all other image parameters
immersive. If the parameter a* were immersive too, the local cyclic
cover would also give length2 there, a contradiction. Its image
y*=pi(a*) is therefore the unique nonimmersion value, although B can
have other multibranch values. At a multibranch value the normalization
parameters stay distinct; only one can be a*.

## 4. A pointwise ledger at every finite target value

Let c(y),b(y),h(y) be the numbers of normalization branches of C,B,H
over y. Count these per prime when H coincides with C or B. Put
delta(y)=1 at y* and zero elsewhere, only in row B. Then

| Row | Omitted ambient lengths over y | Retained source length |
|---|---|---|
| A | c(y) copies2, b(y) copies2, h(y) copies1 | 6-2*c(y)-2*b(y)-h(y) |
| B | c(y) copies2, b(y) copies2 with one increased to3 at y* | 6-2*c(y)-2*b(y)-delta(y) |
| C | c(y) copies2, b(y) copies3 | 6-2*c(y)-3*b(y) |

Every retained point has length1. These formulas include intersections,
self-intersections and exceptional image values. Negative numbers are
forbidden; zero source points are not silently declared impossible.
The equations of C,B,H and the location of y* are still inputs, so the
table is a complete conditional length rule, not a realized global map.

Rows A,C have smooth N everywhere affine. At each omitted ramified
point the completed local map is (s,v)->(s^e,v), so its fibre algebra
is C[s]/(s^e). At a row-B exceptional point the length-three Artin
algebra has only the two possible local algebra types C[s]/(s^3) and
C[s,t]/(s,t)^2. If N is regular there, the fibre is a two-generator
complete intersection and only C[s]/(s^3) can occur. We do not determine
the singularity or algebra type there without further data.

The exact singleton branch-count rules are particularly short:

| Row | Necessary and sufficient integer condition for length one |
|---|---|
| A | (c,b,h)=(0,0,5),(0,1,3),(1,0,3),(0,2,1),(1,1,1),(2,0,1) |
| B | y=y*, with (c,b)=(0,2) or (1,1) |
| C | (c,b)=(1,1) |

These are integer conditions for the displayed fibre formula, not
sufficient geometric existence conditions. In row A, H=C or H=B forces
(c,b,h)=(1,1,1) for a singleton. Thus an off-branch singleton can only
be row A with three distinct target images and five H-branches at zero.

In row B the source cardinality is even at every value other than y*;
at y* it is either3 or1. In row C a singleton is exactly a value where
the two ramified images each have one branch. In particular the branch
counts at y* cannot be substituted by generic e*f alone.

## 5. Inverse and discriminant audit for the residual rows

For a source-finite primitive inverse R(T), let H_C,H_B,H_H define the
images, and let the respective positive source-coordinate pole orders
be p_C,p_B,p_H. With f=1 now proved,

    a_6 = constant*H_C^p_C*H_B^p_B*(H_H^p_H if row A).

If H equals a ramified image, combine its pole exponents. The normalized
trace discriminant is supported with multiplicities (1,1) on C,B in
rows A,B, and (1,2) in row C. PF2OMIT4 therefore gives

    rows A,B: Disc(R)=constant*H_C*H_B*J^2,
    row C:    Disc(R)=constant*H_C*(H_B*J)^2.

The odd-exponent part misses the index-three ramified image B in row C.
All leading pole charges and missing lengths remain separate fields.

At an off-branch singleton of row A, all five omitted points lie on the
one unramified prime D, with its same pole order p=p_H. The sharper
PF2OMIT4 jet bounds are

    ord_0 a_i >= (i-1)*p (2<=i<=6),
    ord_0 a_6=5*p, ord_0 Disc(R)>=20*p, ord_0 J>=10*p.

No actual p, inverse equation or target image polynomial is supplied.

## Replay scope

[verify_residue_ledger.py](verify_residue_ledger.py) checks local fibre
algebras, residue-ramification controls, the forced f1 rows, every bounded
singleton branch pattern, parity, discriminant weights and equal-pole
jet controls. General flat conservation, the image-immersion theorem and
the inherited global r7 obstruction are written arguments/imports.
Local finite controls are explicitly not global Keller realizations.
The full research goal and unrestricted singleton-zero properness remain
open until the actual missing image equations and exceptional values,
or a global exclusion of the remaining complete ledger, are supplied.
