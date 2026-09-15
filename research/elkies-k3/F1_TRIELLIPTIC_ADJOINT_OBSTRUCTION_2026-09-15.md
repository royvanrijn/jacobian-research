# The F1 adjoint system excludes the proposed trielliptic supply

**Written geometric theorem with checked algebraic premises.** A smooth
anticanonical curve R in F1 has no degree-three map to an elliptic curve,
even over an algebraic closure of characteristic zero. Applied to the
halving curves of a 24-I1 elliptic K3, this closes the proposed moving
elliptic supply in the arithmetic-genus-four bisection layer.

Consequently, over a fixed number field there are only finitely many
geometrically integral bisections of arithmetic genus four and normalization
genus one, modulo inherited Mordell--Weil translations. This is finiteness,
not emptiness or an effective enumeration. The
[cubic-contact constructor](CUBIC_CONTACT_AND_TRIELLIPTIC_GATE_2026-09-15.md)
remains valid for individual exceptional contacts. No correlated MW17 cover
is constructed, and larger image genera remain outside this theorem.

## 1. The fixed divisor forced by F1

Write C0 for the negative section and F for the ruling class on F1. Then

```
C0²=-1, C0.F=1, F²=0,
K_Y=-2C0-3F,       R=4C0+6F,       g(R)=9.
```

Let f:R -> P1 be the degree-four ruling and write L=f*O(1). Adjunction gives

```
K_R-2L = (2C0+F)|R.                              (1)
```

All sections of 2C0+F vanish on C0 because its intersection with C0 is -1.
The residual system |C0+F| is basepoint-free. These are all the sections
of the restricted line bundle: in the restriction sequence the kernel is
O_Y(K_Y-2F), whose H0 vanishes and whose H1 is dual to H1(Y,2F)=0.
The latter vanishing follows from the ruling to P1.

It follows that the complete system |K_R-2L| has a fixed divisor C0.R of
length two. This includes tangency: the fixed divisor need not consist
of distinct points. Its space of sections has dimension three.
This restriction argument is the same one used in the
[unique-pencil proof](Q80_HALVING_PENCIL_GEOMETRY_2026-09-15.md).

## 2. A degree-three elliptic map would give a smooth product model

Suppose g:R -> E had degree three, with E elliptic. The maps f and g
generate the function field: the degree of a common intermediate extension
would divide both four and three. Therefore (f,g) is birational onto its
image D in P1 x E.

The [Castelnuovo--Severi bound](https://link.springer.com/article/10.1007/s40993-024-00543-4)
is nine, equal to g(R). More explicitly, D has line bundle

```
O(D)=O_P1(3) external-product A,       deg(A)=4.
```

There is no mixed Picard term because Pic^0(P1)=0. On the product D²=24
and D.K=-8, so its arithmetic genus is nine. Since its normalization R
already has genus nine, D is smooth and R is isomorphic to D. Adjunction
on this model identifies the same intrinsic adjoint bundle as

```
K_R-2L = (O_P1(-1) external-product A)|D.         (2)
```

The following explicit sections show that (2) is basepoint-free. That will
contradict the fixed divisor in (1).

## 3. Three global sections with no common zero

Locally trivialize A on E and use an affine coordinate t on P1. Write D as

```
F(t,e)=a0(e)+a1(e)*t+a2(e)*t²+a3(e)*t³=0.
```

The coefficients are local expressions of global sections of A. On the
finite t chart take

```
q1=a1+a2*t+a3*t²,
q2=a2+a3*t,
q3=a3.
```

For O_P1(-1), passing to the chart s=1/t multiplies a section's local
expression by t. The equation F=0 gives regular expressions there:

```
t*q1 = -a0,
t*q2 = -a1-a0*s,
t*q3 = -a2-a1*s-a0*s².                          (3)
```

These identities hold on D and prove that the q_i define global sections
of (2). They also transform correctly under a change of trivialization of A.

At a finite point, a common zero of q3,q2,q1 forces a3=a2=a1=0; the curve
equation then forces a0=0. At infinity the curve equation is a3=0, and
(3) forces a0=a1=a2=0 if all three sections vanish. But simultaneous
vanishing of all four coefficients at e would make the entire fibre
P1 x {e} a component of D. This contradicts irreducibility and its degrees.
Thus the three global sections have no common zero.

The same bundle cannot be both basepoint-free and have a nonempty fixed
divisor. This contradiction excludes every degree-three elliptic map from R.
The proof does not assume the map is cyclic or defined over the original
number field. It makes no simplicity claim about Jac(R), and does not
exclude elliptic maps of larger degree.

## 4. Consequence for cubic contacts and bisections

For a bisection of arithmetic genus four, the trace-height identity is
h(T)=4*(B.O)+2. It implies geometric primitivity modulo twice MW. The
halving quotient is F1 and its branch curve R is of the type just treated,
as proved in the [cubic-contact note](CUBIC_CONTACT_AND_TRIELLIPTIC_GATE_2026-09-15.md).

The primitive degree-four ruling and genus nine already exclude maps of
degree at most two to genus zero or one, giving finiteness of points of
degree at most two. If infinitely many cubic points existed over a fixed
number field, the [Kadets--Vogt theorem](https://arxiv.org/pdf/2208.01067)
used in that note would force a cubic map to genus zero or an elliptic
curve of positive rank. The genus-zero case has Castelnuovo--Severi bound
six; the elliptic case is excluded above. Hence cubic points are also finite.

There are therefore finitely many rational effective contact divisors of
degree three, including repeated lower-degree points. Contact rigidity
allows at most one integral candidate per divisor. Finally there are
finitely many trace classes modulo twice the inherited MW group, and
translation preserves arithmetic genus, normalization genus and the map
to the original base. This proves the finiteness statement modulo translation.
It also proves finiteness of the literal quadratic extensions realized by
this layer; these are the actual inverse images, not arbitrary constant
twists of their branch supports.

The same cubic-point finiteness applies to the F1 halving curve of the
retained Q80 smooth rational-bisection trace. Its usable degree-two and
cubic branch seeds form finite sets. Degree-four seeds and elliptic maps
of degree four are not excluded by this argument.

## Evidence and next boundary

The [checker](scripts/verify_f1_trielliptic_obstruction.py) verifies the three
chart-gluing identities, invertible coefficient matrices for the two
common-zero arguments, and the intersection and genus numbers. The
[certificate](../artifacts/generated-results/elkies-k3-f1-trielliptic-obstruction-v1/result.json)
distinguishes these finite checks from the written restriction, adjunction,
and arithmetic arguments. No formal verification or external review is
claimed, and no Jacobian computation or point enumeration was needed.

```
.venv/bin/python research/elkies-k3/scripts/verify_f1_trielliptic_obstruction.py
```

Do not search for a degree-three elliptic map on these F1 halving curves:
the proposed supply is impossible. Individual exceptional cubic contacts
are still possible and must pass the retained interpolation test. None is
provided here. The [MW17 two-gain objective](CORRELATED_QUADRATIC_GAINS_2026-09-12.md)
remains open; finite contact sets are not a construction or an exclusion of
all low-genus quadratic covers.
