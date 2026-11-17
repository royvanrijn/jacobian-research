# Small output combinations exclude47 dense low-degree profiles

GLDC5 is a support-independent corollary of the degree-controlled plane
reduction already used in KPAIR1. Among the74 affine-invariant degree flags
strictly below(7,6,4), every one of the47 flags with smallest entry at most2
is impossible for a Keller collision. The27 flags with smallest entry at
least3 remain unclosed by this argument. This is not a full dense census,
a new counterexample, or a global degree/support/dimension optimum.

The subsequent [GLDC6 classical cubic import](CLASSICAL_CUBIC_IMPORT.md)
closes the additional profile(3,3,3). The combined dense status at that stage was
48 excluded/26 remaining. The47/27 partition below is this corollary's own
scope and the baseline for its conditional application counts.

Current later progress is recorded in [GLDC7/8](../quartic-three-audit/README.md):
54 excluded/20 remaining and3D degree[5,7]. The theorem below keeps its own47-profile scope.

## Statement, including the degree range

Let F be a Keller map of affine3-space over a characteristic-zero field,
and put D=max_i deg F_i. Suppose a nonzero constant output combination
R=v^T F has degree s in{1,2}. Then there is no collision if

```
s=1 and D<=100,       or       s=2 and D<=50.
```

The conclusion also holds over the algebraic closure. The theorem imports
only Moh's established plane theorem for both coordinate degrees at most100,
with the same trust boundary as KPAIR1. Moh's theorem is not reproved by
our checker. The numerical caps100/50 are sufficient bounds, not sharpness
claims. The degree-seven census is strictly inside this range.

## A noncritical quadratic has a constant nonzero directional derivative

A nonzero row combination of an invertible Jacobian cannot vanish at any
point. Hence grad R is nowhere zero, and R cannot be constant. Write

```
R(X)=c+ell^T X+(1/2) X^T A X,     A=A^T.
```

If ell belonged to im A, the linear system A*X=-ell would have a solution,
a critical point of R. Thus ell is not in im A. The annihilator of im A
contains a functional w^T nonzero on ell. Symmetry gives A*w=0 and
ell^T*w!=0. Consequently

```
D_w R = ell^T*w != 0
```

is constant. This linear-algebra proof works over the coefficient field,
without a generic-rank assumption, chosen invertible Hessian minor or
omitted singular stratum. A nonsingular A is impossible. For s=1 the same
statement follows immediately from ell!=0.

After a linear source change taking one coordinate direction along w and
a constant target change choosing a rescaled R first, the first output
has the form t+H(u,v), with deg H<=s. Its coefficient of t is1 and all
higher t powers and mixed t terms vanish by the constant derivative.
An affine translation removes any irrelevant constant. Every remaining
output still has degree at most D; invertible affine changes preserve
this maximum degree. The triangular source substitution

```
t -> t-H(u,v)
```

straightens the first output to t. It is a polynomial automorphism and
has determinant1. The resulting map is (t,P(t,u,v),Q(t,u,v)) with nonzero
constant Jacobian. For each fixed t, the pair(P,Q) is a plane Keller map.
Every substituted monomial has u,v-degree at most D*max(1,deg H)<=D*s.
Thus D*s<=100 makes every fibre pair injective by Moh. Colliding source
points would have the same t and contradict this fibre injectivity.
There is no assumed polynomial inverse for the entire parameter family.

For arbitrary characteristic-zero coefficients, all map, transformation
and hypothetical collision coefficients lie in a finitely generated field
that embeds in C. Apply the complex statement after that embedding.
This proves the stated collision exclusion over the algebraic closure.
The written argument, linear algebra and external plane theorem are the
universal proof. Symbolic controls below are not a formalization of them.

## Why the invariant flag is sufficient

For GLDC1's filtered three-dimensional output span V_F, its smallest flag
entry d3 is the least degree of a nonzero element of V_F. An adapted target
basis realizes this degree; no fixed displayed row is required. A constant
nonzero element of V_F would have zero gradient and is impossible for a
Keller map. Thus d3<=2 supplies exactly a linear or quadratic R as above.
The maximum flag entry equals max_i deg F_i in every target basis, because
all basis vectors belong to the same filtered space. Invertible affine
source changes also preserve it. Therefore every flag below(7,6,4) has
D<=7, and the small-component theorem excludes every one with d3<=2.

The checker independently enumerates ordered triples and compares them
with GLDC1's existing degree_profiles_below compiler, without reconstructing
missing historical census artifacts. It finds47 excluded flags:27 with
d3=1 and20 with d3=2. The complete27 remaining flags are explicitly saved
in validation.json. Every one has d3>=3. The original dense quadratic
exclusions are included in these47; this result adds43 full profiles to
the previous explicit four-profile dense closure.

A source/output frame adapted to a collision need not be a degree-adapted
frame. That does not matter: the criterion concerns the whole output span,
and the pairing is affine-left-right invariant. Counting only displayed
component degrees or requiring a displayed quadratic row would miss cases.

## Verification and unchanged boundaries

```
.venv/bin/python research/jacobian/small-component-gate/verify.py
```

The checker pins the KPAIR1 source and its existing validation, confirms the
primary plane source's retained range, checks exact symbolic straightening
and determinant1 for each quadratic Hessian rank0/1/2 normal form, and checks
the unique critical point in rank3. It verifies the monomial substitution
degree bound at the stated100/50 caps, enumerates all74 profiles through
independent code, and writes the47/27 partition plus exact conditional
literature deltas. It does not rerun a sparse census, Groebner basis or Moh's
proof, and does not label these controls an independent universal replay.

The known three-variable degree-seven counterexample has flag(7,6,4),
whose smallest entry4 is outside this criterion. Normalized raw degree at
most7 still has nonlinear-support interval[10,48]; the raw13-term example
is a separate coordinate metric. Ordinary cubic12/36 and cubic-linear89
are unchanged. The global minimum ordinary degree across varying dimensions
remains3, and global dimension/support minima remain open.

After this closure, the classical cubic theorem would add only the flag
(3,3,3), leaving26. The unimported quartic theorem would add four flags
(3,3,3),(4,3,3),(4,4,3),(4,4,4), leaving23. The earlier literature dictionary
counted16 additional flags relative only to the four dense quadratic cases;
that historical reference must not be reused as the current delta.

The separate GLDC6 import now closes333 as stated above. Next attack one
of the26 remaining profiles or a construction that changes the current
positive bounds. Wang-Yan still requires its separate proof audit; its
current additional profiles are433,443,444. Do not repeat a support search
or dense coefficient calculation with d3<=2 in this range.

## Sources

- [KPAIR1's exact reduction and import boundary](../sparse-eight/COORDINATE_GATES.md).
- [KPAIR1 retained exact validation](../sparse-eight/coordinate-gates-validation.json).
- [Moh's primary author account](https://www.math.purdue.edu/~ttm/jacobian),
  history paragraph stating both plane component degrees at most100;
  T.T.Moh, J.reine angew.Math.340(1983),140-212.
- [GLDC1's invariant flag definition](../../foundations/extended-geometry/GLOBAL_LOW_DEGREE_SUPPORT_CENSUS.md).
- [Independent application/controls receipt](validation.json).
