# An affine-source triple-root Cox map with four dicritical components

Orienting the determinant created by minimal pole clearing repairs the
triple-root affine completion.  The result is a polynomial morphism from
affine three-space to a discriminant-type Cox hypersurface.  On smooth
principal charts its hypersurface-residue Jacobian is the constant `1/6`.
The map is birational, one branch component extends finitely, and four
distinct branch components are dicritical.

This is the first construction in the programme with an affine-space source
completion and several dicritical components having distinct target images.
The target is not affine space and is singular where branch components
intersect.

The same ledger has a smaller presentation.  The
[primitive-multiplicity model](PRIMITIVE_MULTIPLICITY_TRIPLE_ROOT_MAP.md)
uses the character `u=ac` and the target
`u^2v=xy(x+1)(x-y+1)`.  It has residue Jacobian `-1`, three distinct
dicritical components, and is unique among monomial factor-splitting
presentations.

Work in characteristic zero.

## 1. Cleared affine coordinates

Use the coordinates `(a,b,c)` and polynomial outputs from the
[triple-root completion](TRIPLE_ROOT_AFFINE_COMPLETION_OBSTRUCTION.md):

\[
\begin{aligned}
A&=-a^2c+3bc+2,\\
C&=-b(bc+1)(a^2c-bc-1),\\
E&=1-a^2c.
\end{aligned}                                       \tag{1}
\]

Their ordinary Jacobian is

\[
 J=-6abc(bc+1)(a^2c-bc-1).                          \tag{2}
\]

Introduce the target orientation coordinate

\[
 R=J.                                                \tag{3}
\]

Set

\[
 x=\frac{A-E-1}{3},\qquad y=1-E.                    \tag{4}
\]

On the source these are simply

\[
 x=bc,\qquad y=a^2c.                                \tag{5}
\]

## 2. Oriented target

Exact substitution gives

\[
\boxed{
 3R^2=
 4C(1-E)(A-E-1)(A-E+2)(A+2E-1).
}                                                     \tag{6}
\]

Let `T_tr` be this hypersurface in `A^4_(A,C,E,R)`.  Equations (1)--(3)
define a polynomial morphism

\[
\boxed{
 \Psi_{\rm tr}:\mathbb A^3_{a,b,c}\longrightarrow T_{\rm tr}.
}                                                     \tag{7}
\]

The branch divisor `R=0` has five distinct components:

\[
\begin{array}{c|c}
\text{component}&\text{short equation}\\ \hline
C=0&C=0\\
1-E=0&y=0\\
A-E-1=0&x=0\\
A-E+2=0&x+1=0\\
A+2E-1=0&E+x=0.
\end{array}                                          \tag{8}
\]

The target is smooth at the generic point of every component and singular
only along loci where branch components meet or the branch arrangement has
higher degeneracy.  Choose a function `h` vanishing on the singular locus
but not identically on any component in (8).  Then `D(h)` is a smooth
affine principal chart retaining the generic point of all five components.

## 3. Constant residue Jacobian

For the defining equation

\[
 H=3R^2-
 4C(1-E)(A-E-1)(A-E+2)(A+2E-1),
\]

one has

\[
 \frac{\partial H}{\partial R}=6R.
\]

On a smooth target chart, use the hypersurface residue form

\[
 \Omega_T=\frac{dA\wedge dC\wedge dE}{6R}.           \tag{9}
\]

By (2)--(3),

\[
 \Psi_{\rm tr}^*(dA\wedge dC\wedge dE)
 =R\,da\wedge db\wedge dc.
\]

Therefore

\[
\boxed{
 \Psi_{\rm tr}^*\Omega_T
 =\frac16\,da\wedge db\wedge dc.
}                                                     \tag{10}
\]

The equality extends across the smooth generic branch points.  Hence the
restriction of (7) over `D(h)` is etale with constant residue Jacobian
`1/6`.

## 4. Birational inverse

On the dense target open

\[
 Cx(x+1)(E+x)\ne0,
\]

the inverse is

\[
\boxed{
\begin{aligned}
a&=\frac{R}{6x(x+1)(E+x)},\\
b&=\frac{C}{(x+1)(E+x)},\\
c&=\frac{x(x+1)(E+x)}{C}.
\end{aligned}
}                                                     \tag{11}
\]

Thus (7) is birational.  Formula (11) also makes the boundary behavior
completely explicit.

## 5. Four dicritical components

At the generic point of any component in (8), `R` is a uniformizer and the
corresponding branch factor has order two.  Applying (11) gives:

\[
\begin{array}{c|ccc|c}
\text{target component}
 &v(a)&v(b)&v(c)&\text{behavior}\\ \hline
y=0&1&0&0&\text{finite: }a=0\\
C=0&1&2&-2&\text{dicritical}\\
x=0&-1&0&2&\text{dicritical}\\
x+1=0&-1&-2&2&\text{dicritical}\\
E+x=0&-1&-2&2&\text{dicritical}.
\end{array}                                          \tag{12}
\]

The four dicritical valuations have four distinct target images.  They are
therefore independent even before adding auxiliary markings.  The fifth
branch extends to the affine source divisor `a=0`.

### Theorem 5.1

The oriented triple-root map (7), restricted to a smooth principal target
chart retaining the generic branch components, is a birational
constant-residue-Jacobian Cox map with four distinct dicritical divisors.
Its polynomial source completion is `A^3`.

## 6. Exact finite-field fibers

The inverse formulas (11) are defined over every field of characteristic
different from `2` and `3`.  Consequently the corresponding dense source
and target opens are isomorphic over every such field and have bijective
sets of rational points over every finite extension.

The complete boundary count is also explicit.  Let `F_q` have
characteristic greater than three.  For each fixed `(A,E,R)`, equation (6)
is linear in `C` away from the zero of its coefficient; when that
coefficient vanishes it forces `R=0` and leaves `C` arbitrary.  In either
case the total count is

\[
\boxed{\#T_{\rm tr}(\mathbb F_q)=q^3.}                \tag{13}
\]

Thus source and target have the same number of rational points.  They are
nevertheless far from bijective.

Using `x=bc`, `y=a^2c`, and

\[
 C=b(x+1)(1-y+x),\qquad
 R=6ax(x+1)(1-y+x),
\]

one obtains the exact target-fiber profile:

\[
\begin{array}{c|c}
\text{fiber size}&\text{number of target points}\\ \hline
0&(4q-5)(q-1),\\
1&q^3-4q^2+5q-2,\\
q-1&3(q-1),\\
q&q-1,\\
2q-1&1.
\end{array}                                          \tag{14}
\]

For example, the unique `(2q-1)` fiber lies over
`x=y=C=R=0`.  The `q`-point fibers have
`x=y=R=0` and `C!=0`.  The three families of `(q-1)` fibers come from
`x=0,C=0`, `x+1=0`, and `1-y+x=0`, with overlaps treated by the displayed
exceptional fibers.

The weighted sum of (14) is `q^3`, as required by the source count, and the
unweighted sum is the target count (13).

Hence:

\[
\boxed{
\Psi_{\rm tr}(\mathbb F_q)
\text{ is not a permutation for any }\operatorname{char}\mathbb F_q>3.
}                                                     \tag{15}
\]

The construction removes the symmetric three-sheet excess on its
birational open, but boundary collisions create a new, exactly quantified
excess.  Permutation behavior still requires a further boundary
modification.

## 7. Reproduction

Run

```bash
.venv/bin/python scripts/verify_affine_source_triple_root_cox_map.py
```

The checker verifies the target equation, residue cancellation, birational
inverse, five branch components, all five valuation vectors, the target
point count `q^3`, and the complete fiber profile (14).
