# An indecomposable Keller map of geometric degree twelve

Work over a characteristic-zero field.  In source coordinates `x,y,z`, put

\[
 u=1+xy,\qquad
 \gamma=1-\frac{12}{11}xy+x^2z.
\]

Define

\[
\boxed{
 F_{12}(x,y,z)=
 \left(
 \frac{10u+u^2-11u^{12}\gamma^{10}}{10x^2},
 \frac{5+u-6u^{11}\gamma^{10}}{5x},
 x\gamma
 \right).
}                                                       \tag{1}
\]

The apparent quotients in (1) are polynomial: the first numerator is
divisible by `x^2` and the second by `x`.  Thus (1) is a concrete element of
`\mathbb Q[x,y,z]^3`, not merely a rational presentation.

> **Theorem.** The map `F_12` has Jacobian determinant one and geometric
> degree twelve.  Its geometric monodromy is `S_12`.  It has no
> decomposition
> \[
> F_{12}=F_2\circ F_1
> \]
> in which both `F_i:\mathbb A^3\to\mathbb A^3` are noninvertible polynomial
> maps.  The same remains true after extending the characteristic-zero base
> field.

## 1. Weighted-seed certificate

Take the sparse admissible seed

\[
 H(W)=\frac{W^2-W^{12}}{10},\qquad c=b_0=1.
                                                               \tag{2}
\]

It satisfies

\[
 H(0)=H'(0)=H(1)=0,\qquad H'(1)=-1,
\]

and

\[
 \kappa=\frac{H''(1)}c=-13,\qquad
 a_0=-\frac{1+\kappa}{2+\kappa}=-\frac{12}{11}.
                                                               \tag{3}
\]

For the weighted construction,

\[
 q(W)=\frac{WH'(W)-H(W)}c
     =\frac{W^2-11W^{12}}{10},\qquad W=u\gamma.
\]

Substitution in the universal weighted formula

\[
 A=\frac{u+q(W)/\gamma^2}{x^2},\qquad
 B=\frac{c+H'(W)/\gamma}{x},\qquad
 C=x\gamma
\]

gives exactly (1).  Weighted polynomiality proves the two divisibilities in
(1), and the weighted-suspension determinant identity gives

\[
 \det DF_{12}=b_0c=1.                                  \tag{4}
\]

For a target `(A,B,C)`, its inverse equation is

\[
 E_{A,B,C}(W)
 =\frac{W^2-W^{12}}{10}-BCW+AC^2=0.                    \tag{5}
\]

On `C\ne0`, every simple root of (5) reconstructs a unique source point.
Consequently the function-field degree of `F_12` is the degree in `W` of
(5), namely twelve.

The change of target parameters

\[
 s=BC,\qquad t=AC^2
\]

is birational on `C\ne0`.  Hence the generic inverse extension is the
incidence extension defined by

\[
 H(W)-sW+t=0.
\]

The [universal symmetric-monodromy theorem](UNIVERSAL_SYMMETRIC_MONODROMY.md)
therefore gives

\[
 \operatorname{Mon}_{\mathrm{geom}}(F_{12})=S_{12}.     \tag{6}
\]

## 2. Atomicity certificate

The [primitive-monodromy atomicity theorem](PRIMITIVE_MONODROMY_ATOMICITY.md)
proves that every characteristic-zero Keller map with primitive geometric
monodromy is absolutely and stably atomic.  Its proof has three steps:

1. primitivity is equivalent to maximality of the point stabilizer, hence to
   absence of a proper intermediate function field;
2. the chain rule forces both polynomial factors of a Keller map to be
   Keller;
3. a geometric-degree-one Keller self-map is a polynomial automorphism.

The natural `S_12` action in (6) is primitive because its point stabilizer
`S_11` is maximal.  Therefore a decomposition of (1) always has an
invertible factor.  In particular, `F_12` is compositionally indecomposable
although its geometric degree is the composite integer twelve.  Moreover,
every stabilization

\[
 F_{12}\times\operatorname{id}_{\mathbb A^r}
\]

has the same property.

## 3. Generator and exact check

The repository generator expands (1) into ordinary polynomials:

```python
import sympy as sp
from jcsearch.weighted import WeightedSeedModel, w

H = (w**2 - w**12) / 10
model = WeightedSeedModel(sp.diff(H, w), c=1, b=1)
F12 = tuple(sp.expand(component) for component in model.mapping())
print(*F12, sep="\n")
```

Run

```bash
.venv/bin/python scripts/verify_weighted_seed_theorem.py
```

to check the seed conditions, the polynomial quotients, the exact map,
`\det DF_{12}=1`, the inverse-incidence identity, and inverse degree twelve.
The `S_12` and maximal-stabilizer steps are theorem-level arguments above,
not conclusions inferred from a bounded computation.
