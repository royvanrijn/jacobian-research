# Target-free inverse ADE by projective birth strata

Date: 2026-09-03.

<!-- status-consumer: EC-K3-INVERSE-ADE-PROJECTIVE-BIRTH-STRATA b4a7edb452e6dcc7 -->

## Outcome

The affine variable in the good-prime birth law can be eliminated exactly.
For a fixed parent core, bridge, graph glue, and good odd prime, every possible
new completion root projects to one explicit point of the projective isotropic
quadric.  Old roots still occupy hyperplane sections.  Consequently

```text
rootless neighbour lines
  = isotropic projective quadric
      minus old-root hyperplane sections
      minus projected scaled dual shells.
```

No marked target core, target isometry class, historical neighbour line, or
surviving-root equality is used.  The construction is exact for nonzero glue
cosets as well as for roots of the child core.

The exact theorem and proof are Theorem H0l.2 in
[`RANK_MUTATION_AND_LIFT_THEOREMS.md`](RANK_MUTATION_AND_LIFT_THEOREMS.md).
The checker
[`scripts/certify_inverse_ade_projective_birth_strata.sage`](scripts/certify_inverse_ade_projective_birth_strata.sage)
exhausts every line in three mass-closed ternary genera.  Across 48
state/prime cases it compares 346 predicted physical root sets with
independent materialized-child root enumerations.  All sets agree, and the
quadric complements predict exactly all 192 rootless lines.  A separate
index-two graph-glue control exhausts six more lines; each has 16 born roots
in the nonzero glue coset, and all six predicted root sets again equal the
independently materialized completed-child sets.  The total is 352 exact set
comparisons.

## The elimination

Let `K` be the positive even parent core, let `p` be prime to the core and
bridge discriminants, and let `N_ell` be the neighbour associated with an
isotropic line `ell` in `K/pK`.  The canonical prime-to-`p` discriminant
identification is

```text
iota:A_(N_ell) -> A_K,
iota([v])=p^(-1)[p*v].
```

Fix a graph class `(a,b)`, and a bridge vector `c` in class `b`.  A born child
root has core part `v` satisfying

```text
iota([v])=a,             v^2=2-c^2.
```

Set `z=p*v`.  Then

```text
z in K^dual,
[z]=p*a in A_K,
z^2=p^2*(2-c^2),
red_p(z)=ell != 0.
```

Conversely every `z` satisfying these four conditions gives the born root
`(z/p,c)` in the completion of `N_ell`.  The converse is the useful part:
write `z=p*x+j*y`, where `y` is an adjusted lift of `ell` and `j` is nonzero
modulo `p`.  The norm equation and `y^2=0 mod 2*p^2` force

```text
<x,y>=0 mod p.
```

Hence `x` lies in the dual of the neighbour kernel and `z/p=x+j*y/p` is
exactly one of the affine-layer vectors.  The existential affine variable
has disappeared from the final condition.

Define the finite scaled shell

```text
S_p(a,c) = {z in K^dual :
              [z]=p*a,
              z^2=p^2*(2-c^2),
              red_p(z) != 0}.
```

Its projective reduction `B_p(a,c)=red_p(S_p(a,c))` is the birth locus for
that physical graph/bridge cell.  The complete birth locus is the union of
these finite sets over the graph glue and the finitely many bridge vectors of
norm at most two.

## Direct rootless and ADE predicates

For an old physical root `(x,c)`, survival is supported on the hyperplane
section

```text
H_(x,c)={ell in Q_p^iso : <x,ell>=0}.
```

Therefore the target-free rootless locus is literally

```text
Q_p^iso
  minus union_(old root lines) H_(x,c)
  minus union_((a,b), c) B_p(a,c).
```

The signs may be quotiented in both unions.  This is the requested forbidden-
stratum description.  The old-root strata have codimension one inside the
quadric; the birth strata are arithmetic zero-dimensional strata.

For a nonempty requested ADE type, attach to each line all old roots whose
hyperplanes contain it and all scaled-shell roots whose projective reduction
is that line.  Their products are computed before constructing a child:

```text
old/old:   <x,x'>+<c,c'>,
old/born:  <x,z>/p+<c,c'>,
born/born: <z,z'>/p^2+<c,c'>.
```

The metric classifier of Theorem H0k then gives the exact ADE type.  Thus an
arbitrary ADE locus is a finite union of incidence cells; rootless is the
special case that is a pure complement.

## What this changes in the planner

The old terminal predicate performed an affine closest-vector query after
each proposed line.  The new orientation is:

1. enumerate each required scaled shell once;
2. project it modulo `p` and deduplicate projective points;
3. compile old-root hyperplanes and birth-point hash tables;
4. solve on the quadric using only incidence and table membership;
5. materialize a child only after the requested metric cell is reached.

This removes the marked-target shortcut and makes the no-birth condition an
actual input-independent constraint.  The small-genus checker is a complete
direct solver: it obtains every rootless line from the complement before its
independent child comparison.

## Complexity and foundry boundary

The elimination is not by itself a uniform speed theorem.  In rank `r`, a
shell at norm proportional to `p^2` can have on the order of `p^(r-2)`
representations, the same scale as the number of points on the projective
quadric.  In rank 15, naively expanding the full core shell `z^2=2*p^2` can
therefore be worse than lazy candidate-wise CVP.  The exact implementation
must select between:

- expanded projected shells when their representation count is small;
- orbit-compressed shells when the parent has useful automorphisms;
- lazy affine CVP when the projected support would be too large.

The theorem solves the missing logical no-birth predicate, not its worst-case
compression problem.  It also does not make the 936 bulk foundry pairs ready:
those rows still lack a compatible source marking and the core/bridge/graph
data needed to define either the affine layers or the scaled shells.  Once
that data exists, no marked target core is needed for rootlessness.

## Replay

```bash
sage -python elkies-k3/scripts/certify_inverse_ade_projective_birth_strata.sage
sage -python elkies-k3/scripts/certify_inverse_ade_projective_birth_strata.sage --check
```

The generated certificate is
[`../artifacts/generated-results/elkies-k3-inverse-ade-projective-birth-strata-v1.json`](../artifacts/generated-results/elkies-k3-inverse-ade-projective-birth-strata-v1.json).
