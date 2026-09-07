# Arithmetic rank seventeen for the degree-two alternate-Q80 pencil (2026-09-03)

<!-- status-consumer: EC-K3-R17-ALTERNATE-Q80-ARITHMETIC-RANK17 a304934727bb3f87 -->

## Outcome

Let `X/QQ` be the K3 surface in Elkies's published rank-17 model, and let
`D` be the exact genus-one bisection labelled `norm12-orbit-11952`.  The
primitive hyperbolic plane

```text
U_alt = <D,O+D>
```

is defined over `QQ`, its orthogonal frame is rootless and integrally
isometric to the alternate-Q80 frame, and the Galois action on
`NS(X_Qbar) tensor QQ` is trivial.  The equivariant Shioda--Tate gate therefore
gives

```text
rank E_alt(QQ(u)) = 19 - 2 - 0 = 17.
```

Renaming the new pencil coordinate `u` as `t` gives the requested statement

```text
rank E_alt(QQ(t)) = 17.
```

This is an arithmetic theorem before equation compilation.  It does not use
or recover seventeen individual sections of the alternate fibration.

The subsequent direct compilation is now complete; see
[`R17_NORM12_ORBIT11952_DIRECT_FIBRATION_2026-09-03.md`](R17_NORM12_ORBIT11952_DIRECT_FIBRATION_2026-09-03.md).
The argument here remains the equation-independent arithmetic-rank
certificate.

## The rational Neron--Severi basis

Write the split lattice of the published rootless fibration as

```text
NS(X_Qbar) = <F,O+F> + R17(-1).
```

The exact published-model certificate supplies seventeen `QQ(t)`-sections.
Their height Gram has determinant `948`, and its determinant-one integral
identification with pinned `R17`, together with the independently proved
geometric Picard rank `19` and saturation, gives rational sections
`Q1,...,Q17` whose frame coordinates are the standard pinned basis.  Hence

```text
B_Q = (F,O,Q1,...,Q17)
```

is an integral basis of `NS(X_Qbar)` consisting entirely of divisor classes
defined over `QQ`.  Its change-of-basis determinant from
`(F,O+F,w1,...,w17)` is one.  In particular

```text
NS(X_Qbar) tensor QQ = NS(X_Q) tensor QQ,
```

so the finite Galois image acts trivially on the full rational
Neron--Severi space.

For reference, if `h_i` is the diagonal entry of the pinned `R17` Gram, then

```text
Qi = ((h_i-2)/2)*F + (O+F) + wi.
```

The offsets `(h_i-2)/2` are

```text
8,33,3,10,1,10,13,5,2,5,11,5,10,6,8,4,17.
```

The certificate retains the full published-to-pinned determinant-one matrix,
so every `Qi` is also an explicit integral group-law combination of the
seventeen published rational sections.

## The alternate marking

The cheapest of the ten degree-two alternate-Q80 witnesses has pinned split
coordinates

```text
w = (-1,-1,3,0,0,1,-1,-1,3,1,1,0,2,-1,1,-2,-2),
D = (3,2,w).
```

Since `w.R17.w=12`, direct intersection in
`U + R17(-1)` gives

```text
D^2=0,  D.F=2,  D.O=1.
```

In the literal rational divisor basis `B_Q`, the two columns of the target
marking are

```text
D = (40,-1,-1,-1,3,0,0,1,-1,-1,3,1,1,0,2,-1,1,-2,-2),
O+D = (40,0,-1,-1,3,0,0,1,-1,-1,3,1,1,0,2,-1,1,-2,-2).
```

Their Gram is exactly

```text
[[0,1],
 [1,0]].
```

The complete 19-by-19 splitting has determinant `-1`.  The positive
orthogonal frame has rank `17`, determinant `948`, zero norm-two vectors, and
is integrally isometric to the alternate-Q80 control.  Thus the geometric
fibre-root space is zero.

## Descent of the pencil and zero

The stored bisection certificate does more than give a Galois-fixed numerical
class: it gives a unique regular member over `QQ` whose branch quartic is
irreducible, squarefree, and disjoint from the source discriminant and trace
denominator.  This is an irreducible smooth genus-one curve over `QQ` in class
`D`.  An irreducible square-zero curve on a K3 is nef.  Since `D` is primitive,
the standard K3 elliptic-pencil theorem gives a base-point-free
two-dimensional space `H0(X,O(D))`; because the divisor and line bundle are
defined over `QQ`, the pencil descends to `QQ`.

The old zero `O` is a rational curve over `QQ` and has `D.O=1`.  Therefore its
map to the new base has degree one, and `O` is the `QQ`-defined zero section of
the new pencil.  This verifies descent of the marked plane
`<D,O+D>` geometrically, rather than inferring it from an abstract lattice
isometry.

## Equivariant Shioda--Tate conclusion

Theorem A2 and Corollary A2.1 of
[`RANK_MUTATION_AND_LIFT_THEOREMS.md`](RANK_MUTATION_AND_LIFT_THEOREMS.md)
apply with

```text
dim V^Gamma = 19,
dim R_alt^Gamma = 0.
```

Consequently

```text
rank MW(E_alt/QQ(u)) = dim V^Gamma - 2 - dim R_alt^Gamma = 17.
```

Equivalently, because both the source and target are rootless and the full NS
representation is trivial, the representation-ring edge check returns zero
change from the published arithmetic rank `17` to the alternate arithmetic
rank `17` on every element of the Galois image.

## Exact replay

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_arithmetic_rank_transfer.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_arithmetic_rank_transfer.sage --check
```

The explicit rational-basis Gram, identity Galois matrix, source and target
`U` columns, rootless-frame checks, fixed spaces, quotient actions, and
all-group-element representation-ring check are stored under
`alternate_q80_application` in
[`../artifacts/generated-results/elkies-k3-arithmetic-rank-transfer-controls-v1.json`](../artifacts/generated-results/elkies-k3-arithmetic-rank-transfer-controls-v1.json).

## Proof boundary and literature

The new base function, a Weierstrass equation, and individual generators of
`E_alt(QQ(u))` remain uncompiled.  Those data are unnecessary for the rank
statement, but will be needed for equation-level specialization work and an
equation-side integral Mordell--Weil basis.

The Neron--Severi quotient description of the Mordell--Weil group and the
Shioda--Tate formula are classical; the Galois action on Mordell--Weil
lattices is due to Shioda.  The tailored equivariant quotient gate used here
is documented with sources in
[`LITERATURE_AND_NOVELTY_MAP_2026-09-03.md`](LITERATURE_AND_NOVELTY_MAP_2026-09-03.md).
Elkies's source paper explicitly states that the published surface has
`NS(X)=NS_Q(X)` of rank `19` and supplies the seventeen rational sections:
[arXiv:2608.25406](https://arxiv.org/abs/2608.25406).
