# Smooth Danielewski resolution of the triple-root map

The final singular line in the canonical triple-root chain has a balanced
resolution which preserves polynomiality from the affine source and keeps
one dicritical component.  The target becomes globally smooth, but it is a
Danielewski threefold of class `L^3+L^2`, not affine three-space.

This yields an exact last-step trilemma:

1. an affine-three target fills the dicritical component;
2. the resolution centered on that component is not polynomial from the
   given source; and
3. the balanced polynomial resolution is smooth but not affine space.

## 1. Balanced modification

The penultimate stage of the
[affine-modification chain](TRIPLE_ROOT_AFFINE_MODIFICATION_CHAIN.md) is

\[
T_2=\{cv=x(x+1)(x-a^2c+1)\}.
\]

Its remaining singular and dicritical component lies over

\[
c=0,\qquad x=-1.
\]

Adjoin the balanced ratio

\[
w=\frac{x(x+1)}c.                                   \tag{1}
\]

On the affine source `A^3_(a,b,c)`, where `x=bc`, this is polynomial:

\[
w=b(bc+1).                                          \tag{2}
\]

The modified target is

\[
\boxed{
T_{\rm D}
=
\{cw=x(x+1)\}
\subset\mathbb A^4_{a,c,x,w}.
}                                                     \tag{3}
\]

The map to `T_2` is

\[
v=w(x-a^2c+1).
\]

## 2. Smooth constant-residue map

Define

\[
\boxed{
\Phi_{\rm D}:\mathbb A^3_{a,b,c}\longrightarrow T_{\rm D},
\qquad
(a,b,c)\longmapsto(a,c,bc,b(bc+1)).
}                                                     \tag{4}
\]

Equation (3) is immediate.

For `H=cw-x(x+1)`, the gradient in `(c,w,x)` is

\[
(w,c,-2x-1).
\]

In characteristic zero these three entries cannot vanish on `H=0`: if
`c=w=0`, then `x` is `0` or `-1`, while `2x+1` is nonzero.  Hence

\[
\boxed{T_{\rm D}\text{ is smooth}.}                  \tag{5}
\]

Use the hypersurface residue form

\[
\Omega_{\rm D}
=\frac{da\wedge dc\wedge dx}{c}.                    \tag{6}
\]

Since

\[
\det\frac{\partial(a,c,bc)}{\partial(a,b,c)}
=-c,
\]

one obtains

\[
\boxed{
\Phi_{\rm D}^*\Omega_{\rm D}
=-da\wedge db\wedge dc.
}                                                     \tag{7}
\]

Thus (4) is a globally defined polynomial constant-residue-Jacobian map
from affine three-space to a smooth affine threefold.

## 3. One surviving dicritical component

On `c!=0`, the inverse is

\[
b=\frac xc.
\]

The boundary `c=0` in (3) has two components:

\[
(c,x)=(0,0),\qquad(c,x)=(0,-1).                     \tag{8}
\]

The first is finite: at `c=0,x=0`, equation (2) gives `w=b`, so `b`
extends.

At the generic point of `c=0,x=-1`, the equation
`cw=x(x+1)` makes `x+1` have order one in `c`, while

\[
b=\frac xc
\]

has order `-1`.  This is a dicritical divisor.  Therefore (4) retains
exactly one dicritical target component.

## 4. Why the smooth target is not affine space

The Danielewski surface

\[
D=\{cw=x(x+1)\}
\]

has two simple degenerate fibers over `x=0,-1`.  Stratifying by `x` gives

\[
[D]
=(\mathbb L-2)(\mathbb L-1)+2(2\mathbb L-1)
=\mathbb L^2+\mathbb L.                             \tag{9}
\]

Consequently

\[
\boxed{
[T_{\rm D}]
=\mathbb L[D]
=\mathbb L^3+\mathbb L^2
\ne\mathbb L^3.
}                                                     \tag{10}
\]

Over `C`, the Hodge--Deligne polynomial excludes an isomorphism with
`A^3`, even after ordinary affine stabilization.

Over every finite field,

\[
\#T_{\rm D}(\mathbb F_q)=q^3+q^2,                  \tag{11}
\]

so the full rational-point map cannot be a permutation.

## 5. The last-step trilemma

More generally, adjoining `g(x)/c` pulls back to

\[
\frac{g(bc)}c.
\]

This is polynomial on the source only if `g(0)=0`.  The modification
centered solely on the dicritical component uses

\[
\frac{x+1}{c}=b+\frac1c,
\]

which is not polynomial on `A^3_(a,b,c)`.

The affine-space reconstruction uses `x/c=b`; it is polynomial, but moves
the `x=-1` component to infinity and ends at the identity target.

To be both source-polynomial and centered on `x=-1`, `g` must be divisible
by `x(x+1)`.  The minimal choice is (1), whose target is the non-affine
Danielewski threefold (3).

### Theorem 5.1

Within last-step modifications `g(x)/c`, no choice simultaneously gives:

1. a polynomial lift from the fixed affine source;
2. an affine-three target; and
3. retention of the final dicritical component.

Allowing more simple roots gives a complete and arithmetically stronger
family.  The
[multi-dicritical Danielewski theorem](DANIELEWSKI_MULTI_DICRITICAL_FAMILY.md)
produces one dicritical affine plane for every nonzero root of `P`, and
nonsplit Galois orbits of those roots give bijections on rational points
over selected finite fields.

## 6. Reproduction

Run

```bash
.venv/bin/python scripts/verify_smooth_danielewski_triple_root_map.py
```

The checker verifies smoothness, residue Jacobian, birational inverse,
boundary valuations, motivic and finite-field counts, and the divisibility
trilemma.
