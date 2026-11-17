# Colored divisor-span obstruction

## Status

This note proves a general necessary condition for Cox fills, masks, units,
and affine modifications recorded on a colored boundary. It is an integral
lattice theorem. It neither constructs a divisor when the condition passes
nor proves affine-space recognition.

The implementation is `ColoredDivisorSpanProblem` and
`colored_proportionality_witnesses` in
[`boundary_package_compiler.py`](../scripts/boundary_package_compiler.py).
The exact regression is
[`verify_toroidal_boundary_feasibility.py`](../scripts/verify_toroidal_boundary_feasibility.py).

## 1. Colored divisor lattices

Let \(C_1,\ldots,C_m\) be the geometric colors of a resolved boundary and
let \(g_1,\ldots,g_n\) be proposed rational divisor generators. Write

\[
 A=(\nu_{C_i}(g_j))\in M_{m,n}(\mathbf Z).            \tag{1.1}
\]

A Laurent monomial

\[
 g^x=\prod_{j=1}^n g_j^{x_j},
 \qquad x\in\mathbf Z^n,
\]

has colored order vector \(Ax\). Let
\(\tau=(\tau_1,\ldots,\tau_m)^t\) be a required derivative, conductor, or
mask divisor.

> **Colored divisor-span theorem.**
> The target \(\tau\) is realizable by a Laurent monomial in the declared
> generators only if
> \[
> \tau\in A\mathbf Z^n.                               \tag{1.2}
> \]
> Equivalently, the class of \(\tau\) in
> \(\mathbf Z^m/A\mathbf Z^n\) has order one. If the augmented matrix
> \([A\mid\tau]\) has larger rank than \(A\), the class has infinite order.
> If the ranks agree but the class has finite order greater than one, a
> multiple of the target is realizable but the target itself is not.

This is immediate from additivity of divisorial valuations. The
determinantal-divisor formula computes the exact class order without
choosing Smith transformations:

\[
 \operatorname{ord}([\tau])
 =
 \frac{\Delta_r(A)}{\Delta_r([A\mid\tau])},           \tag{1.3}
\]

when both matrices have rank \(r\). Here \(\Delta_r\) is the gcd of the
\(r\)-minors. If the augmented rank increases, the order is infinite.

Condition (1.2) is weaker than regular-mask feasibility: a regular monomial
also requires \(x_j\geq0\), appropriate pole support, and any declared
Cartier or gluing constraints. Thus failure is conclusive for the declared
generator architecture, while success only hands the model to the
semigroup and nonlinear gates.

## 2. Proportional colored rows

There is a small local certificate for many span failures. Suppose two rows
of \(A\) have the same primitive direction:

\[
 A_i=a p,\qquad A_j=b p,                              \tag{2.1}
\]

where \(p\in\mathbf Z^n\) is primitive and \(a,b\neq0\). Every vector in the
column span of \(A\) then satisfies

\[
 b\,\tau_i=a\,\tau_j.                                \tag{2.2}
\]

> **Proportional-color obstruction.**
> If
> \[
> b\,\tau_i-a\,\tau_j\neq0,                           \tag{2.3}
> \]
> then \(\tau\) is outside even the rational column span of \(A\). Hence no
> integral, Laurent, or nonnegative monomial in the declared generators can
> realize it. Likewise, a zero generator row with nonzero target order is
> immediately obstructing.

This proves more than bounded enumeration: signs and arbitrarily large
exponents do not help.

The statement also prescribes the next construction. If new generator
columns \(B\) are appended, then for every witness (2.3) at least one new
column \(\eta\) must satisfy

\[
 b\,\eta_i-a\,\eta_j\neq0.                            \tag{2.4}
\]

Otherwise the same row relation persists in \([A\mid B]\). A selector or
unit whose orders vanish on the witness colors cannot repair the defect.
One new column may break several witness classes; the theorem does not claim
that the number of witnesses is a lower bound on the number of generators.

## 3. Compiler gate

A `ColoredDivisorSpanProblem` declares:

- the generator-function columns;
- the target colored orders;
- whether failure is obstructing; and
- a certificate that the named generators exhaust the proposed
  architecture.

The compiler reports the generator and augmented ranks, the exact target
class order, integral-span membership, and one proportional-row witness per
primitive row class. An obstruction is emitted only when both the
conclusive flag and the exhaustive-scope certificate are present. This
preserves the distinction between a theorem about one mask architecture and
a theorem about all possible affine completions.

## 4. The corrected \(F_{20}\) boundary

For the corrected Lecacheux root cover, take the four proposed generators

\[
 (\mathrm{mask}_d,\mathrm{mask}_q,\mathrm{mask}_r,
   \mathrm{selector}_{w-1}).                          \tag{4.1}
\]

The first three repeat the pullbacks of \(d,q,r\). The selector comes from
the global rational \(q\)-conductor cover
\(w^2=(y-5)/(y+3)\); its finite colored orders are all zero. The generator
matrix has rank three, while adjoining the derivative target raises the
rank to four. Thus the target class has infinite order.

The compiler returns six representative proportionality witnesses:

| class | colors \(C_i,C_j\) | row scales \(a,b\) | targets \(\tau_i,\tau_j\) | mismatch \(b\tau_i-a\tau_j\) |
|---|---|---:|---:|---:|
| \(d\) | unramified, index four | \(1,4\) | \(0,3\) | \(-3\) |
| \(q\) | crossing, residual | \(1,1\) | \(1,0\) | \(1\) |
| \(r\) | index two, unramified | \(2,1\) | \(1,0\) | \(1\) |
| triple \(E_1\) | index four, simple | \(4,1\) | \(7,0\) | \(7\) |
| triple \(E_2\) | cluster, simple | \(1,1\) | \(3,0\) | \(3\) |
| \(q\)-\(r\) tangent | \(A\)-cluster, \(B\)-cluster | \(2,2\) | \(2,1\) | \(2\) |

Consequently any broader Cox system that principalizes the derivative must
contain new columns satisfying, respectively,

\[
\begin{gathered}
4\eta_{d,\mathrm{unr}}-\eta_{d,4}\neq0,\qquad
\eta_{q,\mathrm{cross}}-\eta_{q,\mathrm{res}}\neq0,\\
\eta_{r,2}-2\eta_{r,\mathrm{unr}}\neq0,\qquad
\eta_{\mathrm{triple}\ E_1,4}
 -4\eta_{\mathrm{triple}\ E_1,\mathrm{simple}}\neq0,\\
\eta_{\mathrm{triple}\ E_2,\mathrm{cluster}}
 -\eta_{\mathrm{triple}\ E_2,\mathrm{simple}}\neq0,\qquad
\eta_{qr,A}-\eta_{qr,B}\neq0.                        \tag{4.2}
\end{gathered}
\]

The old \(q\)-selector violates none of these equalities because it has zero
finite order. Therefore the next viable \(F_{20}\) ansatz must introduce a
genuinely colored Cox divisor. The earlier \(1458\)-assignment nonnegative
search remains a useful regression, but the unbounded lattice theorem is
the actual obstruction.

## 5. Limits

The theorem does not show:

- that a target in \(A\mathbf Z^n\) lies in the nonnegative semigroup;
- that an abstract colored column is Cartier or principal;
- that proposed masks glue across conductors;
- that inverse-adjugate entries become polynomial; or
- that the resulting source and target are affine spaces.

Those are subsequent gates, applied only after the colored divisor span
contains the target.
