# X1092: discriminant gluing excludes pointed base involutions on every MW17 frame

**Scoped theorem.** On the determinant1092 K3 covered by the retained complete
rootless J2 census, no rootless MW17 elliptic fibration admits an
origin-preserving involution with nontrivial action on its base. This holds
geometrically, including the frames for which no equation has been recovered.
It excludes the parent involution required by the
[order-four conic mechanism](ORDER_FOUR_LIFT_CORRELATED_GAIN_MECHANISM_2026-09-15.md)
on this surface. It does not exclude arbitrary correlated gains, automorphisms
changing the fibration, or other Neron–Severi lattices.

## Necessary integral condition

Let F and O be the fibre and zero-section classes of a rootless fibration on
a Picard19 K3. The classes F and O+F span a unimodular U, so

```
NS = U + (-M),   rank(M)=17,
```

where M is the integral Mordell–Weil height lattice. An automorphism preserving
the fibration and its zero section fixes this U pointwise.

A nonidentity symplectic involution has anti-invariant cohomology lattice
E8(-2), contained in NS, and fixes the transcendental lattice pointwise.
See the [Nikulin cohomology description, sections1.2–1.3](https://sarti.pages.math.cnrs.fr/home/Nikulininvolutionsjune06.pdf).
Consequently its action g on M must satisfy

```
g^2=1,   trace(g)=9-8=1,   g acts trivially on A_M=M*/M.
```

The last condition is essential: unimodularity of K3 cohomology identifies the
discriminant groups of NS and its primitive transcendental complement with
opposite forms. The identity action on the latter forces the identity on the
former. The U summand contributes no discriminant group.

For the column-vector Gram convention H, the exact tests are

```
g^T H g = H,    (g-I) H^(-1) is integral.
```

Indeed M* has coordinate lattice H^(-1)Z^17. A rational eigenspace dimension
count cannot replace this integral test.

## Complete finite check

The [retained census](../elliptic-curves/notes/DET1092_ROOTLESS_J2_CENSUS_2026-09-10.md)
provides all nineteen rootless J2 frame types, their exact Grams, and complete
automorphism-group orders. The present packet records generators, closes each
group by exact integer multiplication, checks preservation of its Gram, and
checks that the resulting cardinality equals the retained full group order.
Thus every isometry in every listed group is inspected. We inherit the census's
coverage and group-order computation; this is not a new independent proof of
Niemeier completeness or a classification of J1 orbits.

| Classes | Group order | Traces of involutions, including identity |
| --- | ---: | --- |
| 1,3–9 | 2 | -17,17 |
| 2 | 4 | -17,-5,5,17 |
| 10,16–19 | 4 | -17,-3,3,17 |
| 11,13,15 | 4 | -17,-1,1,17 |
| 12 | 8 | -17,-3,-1,1,3,17 (with multiplicities in the packet) |
| 14 | 8 | -17,-7,-5,-3,3,5,7,17 |

Exactly five involutions have trace1: one each in classes11,13,15 and two in
class12. Every one acts nontrivially on the discriminant group. For each,
already entry(0,0) of (g-I)H^(-1) is nonintegral:

| Class | Witness entries |
| --- | --- |
| 11 | -20/13 |
| 12 | -263/42, -81/26 |
| 13 | -6/7 |
| 15 | -4/3 |

The certificate retains the entire corresponding matrices, not just these
fractions. There is therefore no compatible symplectic involution fixing a
rootless pointed frame on this surface.

Now suppose eta is an origin-preserving involution acting nontrivially on the
base. If eta is symplectic, the preceding obstruction applies directly.
Otherwise compose it with fibre inversion iota. These commute, since eta
preserves the elliptic group law; both act by -1 on the holomorphic two-form.
Their product is a symplectic involution, remains nonidentity because its base
action is nontrivial, and still fixes F and O. This is again impossible.
This proves the stated theorem for either sign on the two-form.

## Reproducibility and scope

The [input packet](../artifacts/generated-results/elkies-k3-x1092-involution-glue-v1/input.json)
binds the existing census by SHA256 and retains the generators. The
[derived certificate](../artifacts/generated-results/elkies-k3-x1092-involution-glue-v1/result.json)
is checked by [the exact replay](scripts/verify_x1092_involution_glue.py):

```
.venv/bin/python research/elkies-k3/scripts/verify_x1092_involution_glue.py
```

The replay requires SymPy1.14.0 and takes less than one second locally. It
performs no new frame enumeration. An initial batched PARI invocation timed
out at25 seconds without accepted mathematical output; its inputs and failure
record are retained under the packet's preflight directory. Separate bounded
calls supplied the small groups. Their exact closure is the accepted evidence;
the timeout is an implementation event, not a mathematical exclusion.

For a different Picard19 rootless frame, the reusable gate is to seek a
trace1 involution acting trivially on its discriminant group before attempting
an equivariant equation construction. Passing this test would still not prove
geometric realization, the rational Q(i) fixed-base condition, an arithmetic
MW17 marking, or a new section. None is supplied here. The
[original correlated-gain objective](CORRELATED_QUADRATIC_GAINS_2026-09-12.md)
remains open.
