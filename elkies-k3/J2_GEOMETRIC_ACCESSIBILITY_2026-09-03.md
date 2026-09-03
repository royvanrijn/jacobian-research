# Geometric accessibility between `J2` frame classes — 2026-09-03

<!-- status-consumer: EC-K3-H3-ROOTLESS-J2-MINIMAL-ACCESSIBILITY 631f50389e0a3283 -->

## Outcome

There is a useful invariant between `J2` classes, but it must be called a
**one-edge incidence distance**, not yet a metric.  For two frame classes
`[W]` and `[W']` on one K3 surface `X`, define

```text
delta_ell([W],[W'])
  = min { F.F' : F,F' are nef Jacobian fibre classes
                   with frame classes [W],[W'] }.
```

The minimum exists, is symmetric, is zero on the diagonal, and satisfies

```text
[W] != [W']  ==>  delta_ell([W],[W']) >= 2.          (1)
```

No triangle inequality is asserted.  The shortest-path closure of these
one-edge weights is the corresponding route metric.

For the two mass-complete determinant-948 rootless classes on the pinned H3
surface, the exact value is

```text
delta_ell(published R17, alternate Q80) = 2.         (2)
```

Moreover the minimum is attained with the same effective zero curve on both
fibrations.  In both directions the zero-section degree is one.  Thus the
optimal relative marking has the full cross-pairing matrix

```text
                  F'       O'+F'
F                 2          3
O+F               3          2
```

with `O'=O`.  In the notation of Theorem H-1 this is

```text
(d,s,t,z)=(2,1,1,-2),
G_A=A^t*J*A-J=[[12,12],[12,12]].
```

The rank-one bridge is a degenerate relative position: the two embedded
hyperbolic planes share their `(-2)` zero.  This is even cheaper than the
rank-two bridge replacements seen in the foundry route ledger.

This changes the interpretation of the alternate frame's equation frontier.
The historically transported representative has old-fibre degree `11,511`
in the published marking, but that is a cost of one **copy**, not of the
`J2` class.  Ten exact degree-two copies of the alternate class already occur
among the 43 norm-twelve genus-one bisections constructed on the published
R17 equation.  The cheapest stored witness is

```text
norm12-orbit-11952,
w=(-1,-1,3,0,0,1,-1,-1,3,1,1,0,2,-1,1,-2,-2),
D=(3,2,w).
```

Its current equation-complexity proxy has coefficient `L1` norm 13, support
9, and group-addition upper bound 12.  The new construction priority is
therefore to compile the two-dimensional pencil `|D|` directly on the
published rootless equation, not to transport the old Q80 corridor through
its million-scale coordinates.

## Why the lower bound is exact

If two nef primitive isotropic rays have intersection zero, the Hodge index
theorem makes them proportional; primitivity then makes the fibre rays equal.
Changing the zero on one fixed Jacobian fibration does not change its `J2`
frame class.

If `F.F'=1`, the two classes span a copy of `U`.  Equivalently, in a splitting

```text
NS(X)=U+W(-1),       F=e,
```

write `F'=(a,1,w)`, so `w^2=2a`.  The fibre-preserving Eichler transvection

```text
T_x(a,b,w)
  = (a+(w,x)+b*x^2/2, b, w+b*x)
```

with `x=-w` sends `F'` to the other standard isotropic generator of `U`.
Thus the complement of `<F,F'>` is isometric to the source frame; applying the
same argument from the `F'` splitting identifies it with the target frame too.
Intersection one cannot join distinct `J2` classes.  Geometrically, the same
conclusion follows because `(pi,pi'):X -> P1 x P1` would have degree one,
making the K3 surface birational to a rational surface.

This proves (1).  For each of the ten displayed alternate witnesses, the
stored exact chord construction gives an irreducible smooth genus-one curve
with class `D`.  Such a curve has square zero and is nef, so `|D|` is an
elliptic pencil.  Direct intersection gives

```text
D.F=2,       D.O=1.
```

Thus the old zero is a section of `|D|`, and splitting
`<D,O+D>` leaves a rootless determinant-948 frame integrally isometric to the
alternate Q80 control.  This attains (2).  Since a zero curve orthogonal to a
new fibre would be a vertical root, zero cost is at least one when the target
frame is rootless; the displayed shared zero attains that lower bound too.

## A finite classifier from a rootless source

Let the source be rootless, so

```text
NS(X)=U+M(-1),       F=e,
```

and `M` is the Mordell--Weil frame.  Fix a proposed old-fibre degree `d>0`.
Every candidate fibre is

```text
D=(a,d,w),       w in M,       w^2=2ad.              (3)
```

Fibrewise translation by the section indexed by `x in M` acts by the Eichler
transvection above and replaces `w` by `w+d*x`.  Consequently candidates
modulo source-section translation are indexed by the finite set `M/dM`.
For a coset `c`, put

```text
mu_d(c)=min { w^2 : w mod dM=c }.
```

The exact all-section nef gate for a genus-one class is

```text
mu_d(c) >= 2d^2,                                      (4)
```

and the minimum old-zero degree in that translation class is

```text
c_O(c)=mu_d(c)/(2d)-d.                                (5)
```

The divisibility-one gate decides whether (3) belongs to a Jacobian `U`.
For every survivor, split off an integral mate, compute `U_D^perp(-1)`,
enumerate its roots, and run exact integral isometry against the target frame
catalogue.
The finite horizontal-wall test of Proposition C2 completes the nef gate; an
explicit irreducible genus-one representative, as in the 43 R17 controls,
proves nefness directly.

Thus a complete fixed-`d` classifier is:

1. enumerate `Aut(M)`-orbits in `M/dM`;
2. compute `mu_d`, retain the congruence, divisibility, and (4) survivors;
3. apply the finite vertical and horizontal wall gates;
4. split the new `U` and classify its frame by integral isometry;
5. retain the full relative matrix `A`, not only `d`;
6. increase `d` until the requested target class occurs.

The first hit certifies `delta_ell`, because (1) starts distinct-frame searches
at `d=2` and every fixed degree is finite.  Equations should be ranked
lexicographically by at least

```text
(d, c_O, shared-zero penalty, pole order,
 resolved-RR ambient, local transforms, coefficient height).
```

The first three coordinates admit intrinsic marked-lattice minimizations.  The
later coordinates belong to an equation chart and can differ dramatically
between copies of one `J2` class.

## What this does and does not explain

The calculation validates the proposed distinction, but sharpens it:

```text
J2 existence
  != cheap historical marking
  != minimum one-edge accessibility
  != cheap equation compilation.
```

For the alternate determinant-948 class, the `J2` minimum is not the
obstruction: it is optimally two and has an optimally shared zero.  The old
coefficient monster came from transporting a badly placed representative.
The remaining cost is the resolved computation of `H0(X,O(D))` and the
Weierstrass conversion for one of the ten low-degree representatives.

For general foundry outputs, `delta_ell` is still only a lower envelope over
copies.  If the foundry must start from a particular equation marking, record
both the intrinsic minimum and the cost of the currently attached copy.  A
route can then be optimized in the graph whose vertices are **marked**
fibrations and whose projected labels are `J2` frame classes.

## Exact replay

```bash
sage -python elkies-k3/scripts/classify_r17_norm12_isotropic_frames.sage
sage -python elkies-k3/scripts/classify_r17_norm12_isotropic_frames.sage --check
```

The exact record is
[`../artifacts/generated-results/elkies-k3-r17-norm12-isotropic-frame-classification-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-isotropic-frame-classification-v1.json).
It classifies all 43 exact norm-twelve records: 33 give the published R17
class and ten give the alternate Q80 class.  It does not yet compile any of
the ten new elliptic pencils or identify their `J1` surface-automorphism
orbits.

## Literature context

The `J2` classification as classification of frame lattices and its separation
from the finer geometric `J1` classification are standard in the
Kneser--Nishiyama literature; see
[Bertin--Garbagnati--Hortsch--Lecacheux--Mase--Salgado--Whitcher](https://arxiv.org/abs/1501.07484)
and [Braun--Kimura--Watari](https://arxiv.org/abs/1312.4421).
The use of fibre intersection as the degree of an elliptic neighbour and the
equation construction from the corresponding linear system are standard in
explicit K3 work; see [Kumar](https://arxiv.org/abs/1105.1715) and
[Elkies--Kumar](https://arxiv.org/abs/1209.3527).
