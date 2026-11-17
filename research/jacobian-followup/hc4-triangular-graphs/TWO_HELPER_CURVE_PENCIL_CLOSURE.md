# Polynomial curve images close every admitted two-helper overlap

15 November 2026. The remaining two-helper image curves form a
single polynomial pencil. Surjectivity from each source affine
space follows by restricting to a line and using an elementary
integrality argument. A generic forgetful map is flat, even
though it may omit finitely many helper points; flatness prevents
hidden pullback components. The actual disjoint source components
then force the helper curves to be disjoint, and ordinary degree
in the helper plane produces the pencil.

Claims HC4MYIR0/1/2 are recorded in
[MATH_STATUS.json](../../MATH_STATUS.json). The two-helper theorem
works in every admitted stable dimension, not only HC4 size.
Together with HC4MYIQ1 it closes the complete admitted HC4
polynomial-direction overlap family. Residual-direction failures,
nonconstant composite content levels, other output geometries,
both plane endpoints and unrestricted HC4/JC2 remain open.

## 1. Polynomial maps onto an affine curve are already surjective

Let `Omega` be algebraically closed and let a morphism
`phi:A^N->C` be dominant, where `C` is an integral affine
curve. Then `phi` is surjective on geometric points, and

```
O(C)*=Omega*.                                      (1)
```

Choose an affine coordinate on `C` whose pullback is nonconstant.
There are two source points where its values differ; restriction
to the affine line joining them is nonconstant. The resulting
map `A1->C` is dominant because its image is not a point.
It identifies `A=O(C)` with a subalgebra of `Omega[t]`.
Some element of `A` is a polynomial `a(t)` of positive degree,
with nonzero constant leading coefficient. Dividing
`a(S)-a(t)` by that coefficient gives a monic equation for
`t` over `A`. Thus

```
Omega[t]=A[t] is finite integral over A.
```

Lying over makes this restricted line map surjective onto `C`.
It lies inside the original image, proving surjectivity of `phi`.
Finally the dominant map injects `O(C)` into a polynomial
source ring. A unit and its inverse map to polynomial units,
so both are scalar. This proves (1).

This argument needs neither smoothness nor normality of `C`.
It does not assert that a dominant map from affine space onto
a higher-dimensional variety is surjective. Primary context is
[Arzhantsev, arXiv:2209.08607v2, Section 1](https://arxiv.org/html/2209.08607v2),
including the scalar-unit observation and Lemma 1(c) on curve
images. The exact assertion used here has been proved directly;
we do not import the curve-normalization classification, flexibility
theorems, Rentschler, Abhyankar--Moh or a cancellation theorem.

## 2. An image inside the multiplier divisor already closes (HC4MYIR0)

Put `k=C(T)`, `F=k(V)`, `L=C(s,t)`, with `T=f,V=g`
for a complex plane Keller pair. Take nonzero `M in k[X]`
and `P in k[X,Z]` with `p=#X>=1,q=#Z>=1`, and assume
the residual polynomial-direction condition of HC4MYIO1:

```
P modulo every geometric factor a_j of M
has positive joint Z degree.                       (2)
```

If a genuine closed stable hypersurface embedding carries the
pair through `(T,M Y-P)` and `M` vanishes identically on its
image, then the pair is invertible. This generic conclusion
works in every admitted stable dimension and does not require
the global content-fibre condition.

Over `Omega=F_bar`, choose the geometric generic image component
corresponding to the original source embedding `L->Omega`.
Its stabilizer is `H=Gal(Omega/L)`. Since `M` vanishes there,
it is a removed prime divisor `D_(j,i)` in the exact HC4MYIO1
row calculation. Its uniquely determined multiplier row `a_j`
is therefore fixed by `H`.

Choose a finite Galois extension `E/k` containing coefficients
of all normalized geometric factors of `M`. Row permutation
factors through `Gal(E/k)`. The source regular-parameter theorem
HC4MYHX0 makes `L` linearly disjoint from `E`, so `H` has
the full target image in this finite parameter group. A row
fixed by `H` is consequently fixed by the full target group
`G=Gal(Omega/F)`.

The subgroup `Gal(Omega/L E)` is transitive on all components
within row `j`, by HC4MYIO1 and its independent coefficient-field
proof. It fixes the selected source component, so that row has
exactly one component. The unique component of a target-fixed
row is target-fixed. Its target orbit has size one, whereas
the actual source components have transitive orbit of size
`[L:F]`. Hence `[L:F]=1`, and the retained birational Keller
case supplies the inverse.

This checks the multiplier-zero image before any localization.
It does not assume the row was originally defined over `k`,
that the source extension is Galois, or that a finite parameter
cover is itself an affine-plane source.

## 3. Every two-helper polynomial-direction overlap closes (HC4MYIR1)

Now take `X=(X1,X2)`, arbitrary `q>=1` root variables `Z`,
and polynomials

```
M(T,X) in C[T,X], M!=0,
P(T,X,Z) in C[T,X,Z].
```

Assume (2), and the exact content condition:

```
P(c,X,Z) is nonconstant closed over C
at every zero c of the coefficient content of M.    (3)
```

Every genuine closed embedding carrying a plane Keller pair through

```
Pi=(T,M(T,X)Y-P(T,X,Z)), r>=q+1,                   (4)
```

forces the polynomial inverse. Include `u=r-q-1` unused variables
`W`. All multiplier, root, coefficient, source, parameter and
embedding degrees, root-variable counts and stable dimensions are
unrestricted. No monicity, Bezout leading certificate, primitive
content, generic positive filtration or separate forgetful-flatness
hypothesis is imposed.

### 3.1 Reduce to an invariant helper-plane equation

Under hypothetical nonautomorphy, HC4MYIP2 gives the actual generic
image equation

```
E_F in F[X1,X2].                                   (5)
```

The same hypotheses give global faithful flatness over `C[T,V]`,
geometrically integral affine height-one fibres, a smooth geometric
generic ambient domain, and scalar geometric generic units. The
actual generic image has geometric components `A^r`, disjoint
and indexed by the `[L:F]` embeddings of `L` into `Omega`.
The equation `E_F` is irreducible in the generic ambient ring.
It is also irreducible in `F[X]`: a nonconstant helper factor
could not become an ambient unit, since those units are `F*`.

Section 2 closes the case where `M` vanishes on the image.
Assume henceforth it does not. Then `E_F` and `M` have no
common geometric helper factor. Indeed any common factor would,
by Galois conjugacy of the factors of irreducible `E_F`, make
`E_F` divide a power of `M`; the generic source quotient is
a domain, so `M` would vanish identically there.

### 3.2 The generic forgetful map is flat; its finite holes are allowed

Write `P=sum_alpha c_alpha(T,X) Z^alpha`, with constant coefficient
`c_0`. Over `k_bar`, the common zero locus

```
Sigma=V(M, c_alpha for |alpha|>0) subset A2_X       (6)
```

has no curve component by (2): an irreducible equation of such
a component would be a multiplier factor on which every positive
root coefficient vanishes. Thus `Sigma` is finite or empty.
Its geometric points have coordinates in `k_bar`. At any such
point `x`, `c_0(x)+V` is nonzero because `V` is transcendental
over `k_bar`. The same set after extension to `Omega` therefore
has empty fibres under the generic forgetful map

```
psi:Spec R_Omega -> A2_X,
R_Omega=Omega[X,Y,Z,W]/(M Y-P-V).                   (7)
```

At a point outside `M=0`, eliminate `Y`: its fibre is `A^(q+u)`.
At a point of `M=0` outside `Sigma`, `P(x,Z)+V` is nonconstant,
so its nonempty principal zero scheme has pure dimension `q-1`.
The free `Y,W` variables give fibre dimension `q+u`. Fibres
may be reducible or nonreduced; the dimension statement is what
is used here. There is no fibre with an extra full root-variable
dimension, since (6) has just been excluded by transcendental `V`.

The source of (7) is smooth of dimension `q+u+2`, and the base
is the regular helper plane. Every nonempty closed fibre has
the expected pure dimension `q+u`. The retained
[miracle-flatness lemma, Stacks 00R4](https://stacks.math.columbia.edu/tag/00R4),
applied to these local rings, makes `psi` flat. Explicitly, over
a closed helper point the two helper-coordinate differences form
a regular sequence in the Cohen--Macaulay local source ring by
the dimension drop two. The dimension equality gives flatness;
height-one base points use torsion-freeness over the DVR. This
checks all local points. The map need not be faithfully flat
on the whole helper plane: the finite set `Sigma` can be omitted.

### 3.3 Helper curves pull back to exactly the actual image components

Factor `E_F` in `Omega[X1,X2]` into distinct irreducible helper
polynomials `g_i`, initially allowing multiplicities, and put
`C_i=V(g_i)`. No `g_i` divides `M`. Base change of the flat
map (7) makes

```
S_i=R_Omega/(g_i)
```

flat over the domain `Omega[C_i]`. At the generic point of
`C_i`, `M` is a unit, so the generic fibre is the polynomial
domain `Omega(C_i)[Z,W]`. Flatness makes `S_i` torsion-free
over `Omega[C_i]`, hence injective into that generic fibre.
Thus `S_i` is a domain: no component supported above a special
helper point, and no hidden nilpotent, has survived localization.

The pulled-back primes are distinct, since their generic images
are distinct curves. The actual geometric generic source quotient
is reduced, so all the multiplicities in `E_F` are one. These
primes are exactly its geometric generic image components.
Consequently each `S_i` is the actual polynomial ring in `r`
source variables. This proves a bijection

```
helper curves C_i <-> source components A^r_i,
number of curves = [L:F].                          (8)
```

### 3.4 Surjectivity and units make the helper curves a pencil

The source component `A^r_i` maps dominantly to `C_i` by the
generic-fibre calculation. Section 1 makes this map surjective,
and makes `Omega[C_i]*=Omega*`. If `C_i` and `C_j` met,
surjectivity would supply a point on the first source component
above their intersection. At that point both `g_i` and `g_j`
would vanish, contradicting disjointness of the actual source
components. Thus the helper curves are pairwise disjoint.

All `g_i` are Galois conjugate up to scalar. Normalize their
common last nonzero monomial coefficient to one. They have equal
ordinary total degree in `X1,X2`. Disjointness and scalar units
in `Omega[C_i]` give

```
g_j-c_ij=b_ij g_i, c_ij in Omega*.
```

Ordinary polynomial degrees add, so `b_ij` is scalar; the common
normalization makes it one. Their nonconstant parts agree.
Galois descent gives a nonconstant `h0 in F[X1,X2]` and labels
`alpha_i` with

```
E_F=Q(h0), Q in F[S], deg Q=[L:F].                 (9)
```

This degree argument is in the polynomial **helper plane**, not
an unproved positive filtration of the generic ambient ring.
The helper coordinates are algebraically independent over `F`
in that ambient ring after inverting `M`, so `h0` is
transcendental over `F` there.

HC4MYIQ0 applies to (9), using the already proved global
faithful flatness and height-one integral fibres. It gives
`deg Q=1`. The bijection (8) now makes `[L:F]=1`, and the
birational Keller theorem gives the polynomial inverse. This
contradicts nonautomorphy and proves (4).

## 4. The admitted HC4 overlap family is closed (HC4MYIR2)

At HC4 size `r=2`, any polynomial-direction overlap with
`p=#X>=1`, `q=#Z>=1` and `u` unused variables has

```
p+q+u=3.
```

The one-helper cases close by HC4MYIQ1. The only other case
has `p=2,q=1,u=0` and closes by HC4MYIR1. Therefore every
closed embedding carrying a plane Keller pair through this
entire HC4-size family, under (2)--(3), forces the inverse.
All degrees, moving factors, repeated multiplier factors,
overlapping coefficients and admitted nonprimitive content
are unrestricted. Ambient and target polynomial automorphism
orbits inherit the conclusion as before; this does not assert
nonlinear source Hessian invariance.

The earlier HC4MYIQ2 invariant plane-curve residual is thus
retained as an intermediate reduction, no longer an open case
inside this admitted family. No assertion is made that arbitrary
HC4 graph embeddings admit observables of this form satisfying
(2)--(3). Residual-direction failures, nonconstant composite
content levels, other output geometries, compatible arbitrary-graph
observables, actual endpoint representatives, both plane endpoints
and unrestricted HC4/JC2 remain open. In larger stable dimensions,
three or more multiplier helpers are a separate geometric question.

## 5. Controls and exact boundaries

The curve-surjectivity statement allows singular curves: for
example `(s,t)->((s+t^2)^2,(s+t^2)^3)` surjects onto the cusp.
It does not require identifying the cusp with `A1` or importing
an embedding theorem.

The dimension-two helper restriction has a real role. With two
helpers, `M=X1`, `P=X2 Z` has the single coefficient-drop point
`(0,0)`; its generic forgetful fibre is empty there because
`V` is nonzero. With three helpers, `M=X1`,
`P=X2 Z+X3` has a coefficient-drop line `X1=X2=0`.
The point on it with `X3=-V` gives a full `Y,Z` fibre, larger
than the generic dimension. The flat-helper-plane argument cannot
be extended to that example by the same finite-point reasoning.
It is a boundary control, not a Keller realization.

Nor does affine-space dominance imply surjectivity onto surfaces.
The polynomial map `(a,b)->(1+ab,a+b^2+ab^3)` has image
`A2 minus {(0,0)}`. If its first coordinate is zero, its second
is `a!=0`; if the first coordinate is nonzero, a cubic equation
in `b` gives a preimage over an algebraically closed field.
This example is also in the inspected primary context. No
surface-image or stable-cancellation inference is used here.

For every `m>=2`, a genuine HC4 two-helper control is

```
M=1+X1 X2,
P=(X1+1)Z^m+T Z+X1 C(T,X2,Z), deg_Z C<m,
j(s,t,u,v)=(s,0,u,t+v^m+s v,v).
```

It carries the identity pair, with source inverse
`(T,Y-Z^m-T Z,X2,Z)`. The nondegenerate ambient quadratic
`T X2+Y Z+X1^2` pulls back to `s u+t v+v^(m+1)+s v^2`.
Its actual four-variable Hessian determinant is one. This tests
overlapping terms, arbitrary root degrees, an actual closed
embedding and actual unit Hessian, without constructing a
nonautomorphic example.

The [source context](sources/two-helper-curve-context.json),
[checker](verify_two_helper_curves.py),
[receipt](receipt-two-helper-curves.json) and
[validation](validation-two-helper-curves.json) keep the written
all-degree proof separate from exact controls and retained byte
pins. The generic helper flatness uses the retained exact
miracle-flatness import. No new unproved theorem, curve or action
classification, stable cancellation, formal verification,
independent replay, external review or novelty priority is claimed.
