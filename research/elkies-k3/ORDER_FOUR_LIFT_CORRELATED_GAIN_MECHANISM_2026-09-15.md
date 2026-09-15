# An order-four lift forces an orthogonal pair from one new direction

**Conditional construction theorem, with an explicit K3 control.** Suppose
an origin-preserving involution eta of an elliptic fibration lifts to an
automorphism alpha over a quadratic cover, and alpha² is the quadratic
deck involution sigma. Then the anti-invariant Mordell--Weil space is a
Q(i)-vector space. In particular its rank is even, and any nontorsion
anti-invariant P gives the independent pair

```
P, alpha(P),       height Gram = diag(h(P),h(P)).
```

Over Q, this can furnish a genus0 cover when the two fixed points of eta's
base involution are conjugate over Q(i). It cannot furnish a genus1
quadratic cover. Thus the condition reduces a two-gain construction to
one new direction on a distinguished conic, provided the required
involution exists on an actual MW17 parent.

The control below realizes the mechanism, but its parent has arithmetic
rank at most10. The four retained MW17 parents have already failed the
[symmetry gate](MW17_J_MONODROMY_AND_SYMMETRY_BARRIER_2026-09-15.md).
The [complete X1092 frame gate](X1092_ROOTLESS_INVOLUTION_GLUE_GATE_2026-09-15.md)
now excludes the required parent involution on all nineteen rootless frames
of that surface. The canonical involution of the determinant622 Inose source fails a
separate fixed-point-field gate, even after changing equivariant fibrations.
The [requested MW17 construction](CORRELATED_QUADRATIC_GAINS_2026-09-12.md)
therefore remains open.

## 1. The representation and height argument

Let E/Q(t) be nonisotrivial, let C -> P1_t be a quadratic cover, and
suppose alpha is an origin-preserving automorphism of the pullback
elliptic fibration with alpha²=sigma. Work modulo torsion. On

```
V^- = {P in E(Q(C)) tensor Q : sigma(P)=-P}
```

one has alpha²=-1. Hence Q[alpha] is Q(i), and dim_Q V^- is even.
For a nonzero P in V^-, a relation aP+b alpha(P)=0 with rational a,b
implies (a²+b²)P=0 after applying a-b alpha, so a=b=0.

The Shioda height pairing is preserved by an origin-preserving fibration
automorphism, since it preserves the intersection calculation and fibre
correction terms. Consequently

```
<P,alpha(P)> = <alpha(P),alpha²(P)> = -<alpha(P),P> = 0.
```

Both diagonal heights equal h(P)>0. Anti-invariant points are orthogonal
to the inherited subgroup. This proves independence modulo that subgroup.
Given any verified new direction Q, one may use P=Q-sigma(Q); the phrase
new direction means nonzero in the quotient after tensoring with Q,
not merely a point outside a chosen finite-index section subgroup.

## 2. The rational conic and its arithmetic condition

Normalize the base involution to phi(t)=-1/t. Take

```
C: w²=t²+1,
beta(t,w)=(-1/t,w/t).
```

Then beta preserves C and beta²(t,w)=(t,-w). The conic is rational:

```
t=(s²-1)/(2s),   w=(s²+1)/(2s),   s=t+w.
```

On s, beta is (s-1)/(s+1), a rational automorphism of order4.
If eta is an origin-preserving involution of E over phi, its pullback
combined with beta gives alpha²=sigma. All maps extend across the omitted
affine points on their smooth projective models.

This base-field condition is necessary for any such rational conic.
An order-four element of PGL2(Q), represented by B, has eigenvalue ratio
+/-i. Thus tr(B)²=2det(B), and the discriminant of its fixed-point equation
is -tr(B)². Its fixed points lie over exactly Q(i). Its square has the
same two fixed points. Their distinct images under the quadratic quotient
remain conjugate over Q(i), so the induced involution on P1_t also has
fixed-point field Q(i). Conversely, a rational involution with this fixed
field is conjugate over Q to -1/t and admits the displayed construction.
This is a condition on the actual rational base involution, not merely
on a geometric order-two lattice isometry.

There is no genus1 version over Q. If beta²=sigma on a genus1 curve with
quotient C/sigma=P1, sigma acts by-1 on the one-dimensional Q-space of
regular differentials. But beta acts on that space by some lambda in Q,
which would require lambda²=-1. The contradiction does not depend on
whether the genus1 curve has a rational point.

## 3. Explicit positive mechanism control

Set

```
A=-t^6+5*t^4-t^2,
B=t^12+t^10+3*t^7-3*t^5+t^2+1,
E: y²=x³+A*x+B.
```

The involution

```
eta(t,x,y)=(-1/t,x/t^4,y/t^6)
```

preserves E and its zero section. On C:w²=t²+1 take

```
P=(t, w*(t^5+1)),
Q=(-t^3, w*(1-t^5)).
```

Exact polynomial identities put both sections on E. Applying the lifted
symmetry to P gives Q; the orbit is P,Q,-P,-Q. The degree24 discriminant
is squarefree modulo13 and coprime to t²+1. The infinity fibre is smooth,
so E is a24-I1 elliptic K3 and the cover branches only over smooth fibres.
The pullback has chi=4 and irreducible fibres. Both sections have no
intersection with zero, including in the minimal infinity charts; their
heights are8. Thus their Gram matrix is diag(8,8).

An independent intersection check gives the same cross term. The finite
abscissa difference is t(t²+1). At t=0 there are two unramified common
points, each of intersection multiplicity1. At each branch t=+/-i,
the ordinate difference 2*w*t^5 has order1 in w, giving one intersection.
The infinity values of the two ordinates are opposite and nonzero, so
there is no intersection there. Hence P.Q=4 and <P,Q>=chi-P.Q=0.

This control fails the high-parent-rank condition. At13, explicit finite
fibre enumeration, also checked by Legendre sums, gives

```
#S(F13)=200,   trace(Frob|H²)=200-1-13²=30.
```

If r rational divisor classes are independent, they supply r eigenvalues
13 in H²; the other22-r eigenvalues have absolute value13. Therefore
30 >= (2r-22)*13, so r<=12. The two trivial fibration classes then give
arithmetic MW rank at most10. This is an unconditional upper bound,
not an exact rank calculation. No coordinate change or extra section on
this parent can upgrade it to the required arithmetic MW17 parent.

## 4. A fixed-point-field filter before lattice or equation work

Suppose eta has a fixed point defined over a number field K that does not
contain i. If an eta-equivariant rational elliptic fibration existed with
the required order-four conic lift, the image of that point would be a
K-rational fixed point of its base involution. Section2 says those fixed
points have field Q(i), a contradiction.

Thus a single such fixed point excludes this lift mechanism for **every**
eta-equivariant rational fibration. This does not exclude other involutions
or arbitrary new-section constructions. It also applies to Q-conjugates
of eta, since a Q-automorphism transports the fixed K-point to another
fixed K-point.

Apply this to the retained
[determinant622 Inose source](DET622_INOSE_RATIONAL_MARKING_SOURCE_2026-09-15.md).
Its untwisted equation has the form

```
y²=x³+a*t^4*x+t^5*(t²+b*t+k),
a=-3IJ, b=-2IJ², k=I²J³.
```

Its canonical symplectic involution is

```
eta(t,x,y)=(k/t, k²*x/t^4, -k³*y/t^6).
```

The zero-section points at t=+/-sqrt(k) are fixed, defined over
K=Q(sqrt(J)). The retained I and J are nonzero and negative. The checker
verifies that -J is strictly between two consecutive integer squares.
Hence K is imaginary quadratic but is not Q(i). The fixed fibres are
smooth as well, by an exact residual-discriminant gcd.

The same involution and fixed zero-section points exist on every constant
quadratic twist of this model. Therefore the still-uncomputed rationalizing
twist does not affect this obstruction. Changing to an equivariant rootless
frame, if one exists, cannot rescue this particular order-four mechanism.
No rootless-frame enumeration or primitive-section reconstruction is needed
for that conclusion. Other involutions are outside this fixed-field argument. Subsequently, the
[independent global root obstruction](DET622_GLOBAL_ROOT_OBSTRUCTION_2026-09-15.md)
proved that every frame of this determinant622 NS lattice is rootful, so this
surface cannot supply a MW17 parent at all. The fixed-field calculation remains
a retained mechanism filter, not an active equivariant-frame search.

For orientation, a symplectic involution preserving a rootless pointed
frame on a Picard19 K3 would split its MW17 space into dimensions9 and8:
the anti-invariant lattice is E8(-2), contained in NS, by the standard
[Nikulin cohomology description](https://sarti.pages.math.cnrs.fr/home/Nikulininvolutionsjune06.pdf),
so it has dimension8,
and the fixed fibre/zero classes consume two of the11 invariant NS
directions. Composing with fibre inversion reverses those MW signs,
giving the eigenspace pattern needed by the proposed8+9 descent.
This alone does not establish a rational-surface quotient. Such a lattice split
alone supplies neither the arithmetic Q(i) condition nor the new section.
The fixed-point-field test must come first for the cyclic-lift route.

## 5. Evidence and status

The [checker](scripts/verify_order_four_gain_mechanism.py) and
[certificate](../artifacts/generated-results/elkies-k3-order-four-gain-mechanism-v1/result.json)
verify the coordinate identities, rational conic parametrization, control
smoothness and point counts, and the exact Inose nonsquare/fixed-fibre
conditions. The computation took about0.05 seconds under20 CPU seconds
and1GiB with SymPy1.14.0. Replay with:

```
.venv/bin/python research/elkies-k3/scripts/verify_order_four_gain_mechanism.py --record /tmp/order-four-gain-replay.json
```

The representation, height, differential and cohomological arguments are
written mathematics. No formal verification, independent whole-proof
replay or external review is claimed. The mechanism is a conditional
sufficient theorem; its explicit example has a low-rank parent. An actual
arithmetic MW17 parent meeting its hypotheses has not been supplied.
