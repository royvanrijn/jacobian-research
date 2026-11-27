# Ordinary sextic inversion for linear one-parameter actions

**Theorem SexticKellerDiagonalTorusInversion.**
Let \(F:k^3\to k^3\) be a characteristic-zero polynomial Keller map
of ordinary degree at most six. Suppose \(F\) is equivariant for
one-parameter algebraic torus actions on source and target, with a
nontrivial source action. Assume the actions are diagonal in affine
coordinates over the algebraic closure of \(k\). Then \(F\) is an
automorphism. The least ordinary degree of a counterexample in this
class is seven.

The proof uses the ordinary total degree in the given affine frames.
The target spectrum equals the source spectrum by differentiation at
the fixed origin, as in [Shaska, Lemmas 2.1–2.2](https://arxiv.org/html/2607.20210v1#S2).
Shaska's Theorem 3.2 treats positive weights without a degree bound.
The quotient calculations of
[Kistner and Shaska, Section 2](https://arxiv.org/html/2608.02863v1#S2)
organize known mixed-sign constructions. The finite reduction here
allows arbitrary integral weights and uses the ordinary sextic cap.

## Normalization and weights

Work over the algebraic closure. Translate the target by \(F(0)\),
a fixed point, and compose on the target with \(JF(0)^{-1}\).
The normalized map has constant determinant one and is the identity
in degree one. Its two weight lists now agree, since \(JF(0)\)
intertwines the linear actions. All transformations so far are affine.

If the weights have one sign, reverse them if necessary to make all
nonzero weights positive. The weight-zero outputs depend only on
the weight-zero variables, of which there are at most two. Differentiating
in order of increasing weight gives a block triangular Jacobian.
Its zero-weight block is the Jacobian of a Keller map of ordinary
degree at most six, which is invertible by the
[planar bound](LOW_DEGREE_COORDINATE_KELLER_MAPS_IN_THREE_VARIABLES.md)
when there are two zero weights. The one-variable case follows
from its constant nonzero derivative. At a positive weight,
each output is linear in the inputs of that weight with
coefficients in the weight-zero variables, plus terms from lower
weights. The determinant of each diagonal block is a unit since their
product is one. Each positive block can therefore be solved
polynomially after the lower blocks.

In the mixed-sign case reverse signs and permute coordinates to use
primitive weights \((a,-b,-c)\) with \(a,c>0,\ 0\le b\le c\).
If the positive-weight output is \(x\), inversion follows from
the [bounded coordinate theorem](LOW_DEGREE_COORDINATE_KELLER_MAPS_IN_THREE_VARIABLES.md).
We can therefore assume it has a nonlinear term. If \(b=0\),
the weights are distinct and the \(y\)-axis is fixed. The determinant
on that axis is the product of three diagonal polynomials in \(y\).
Their values at zero are one, so they are all one. Also \(F_2=y\)
on this axis. Hence the normalized supports omit \(xy^j\) in
the first component, \(zy^j\) in the third, and \(y^j\) in the
second for all corresponding nonlinear terms.

For a nonlinear monomial with exponent \(e\) in component \(i\),
call \(e-e_i\) its resonance vector, where \(e_i\) is the \(i\)-th
standard basis vector. Each vector is perpendicular to the weight
vector. A nonlinear term of the positive-weight output has a
nonnegative resonance vector with at least two positive coordinates.
Indeed, a monomial in that output requires a positive exponent of \(x\).
If it has exponent one in \(x\), every other variable it uses must
have weight zero, and the fixed-axis calculation excluded all those
nonlinear terms.

If the resonance vectors of the actual nonlinear terms span only one
direction, take its primitive integral generator \(v\) in the
nonnegative direction. No negative multiple of \(v\) can be a resonance
vector, since its negative entries in at least two coordinates cannot
be made nonnegative by adding one \(e_i\).
Thus \(F_i=x_i f_i(M)\), with \(M=\prod_i x_i^{v_i}\) and
\(f_i(0)=1\). Its determinant is
\[
 \prod_i f_i(M) + M\sum_i v_i f_i'(M)\prod_{j\ne i} f_j(M).
\]
At its largest possible power of \(M\), the coefficient is the product
of the leading coefficients of the \(f_i\), times
\(1+\sum_i v_i\deg f_i\). This is nonzero in characteristic zero.
Constant determinant therefore forces every \(f_i\) constant,
excluding the assumed nonlinear term.

## Finite checks for a two-dimensional resonance span

Otherwise a positive-output resonance vector and a second,
independent vector determine the weights up to sign and common factor
by their cross product. Both vectors come from monomials of ordinary
degrees two through six, so taking all such pairs yields an exhaustive
finite set. For each candidate weight use every compatible monomial
up to degree six, with the fixed-axis exclusions above when \(b=0\).
Treat the coefficients of degree at least two as parameters. For each
choice of the first nonzero positive-output coefficient, set its
predecessors to zero and that coefficient to one. This normalization
is possible by diagonal conjugation over the algebraic closure:
for the chosen term the exponent of \(x\) in its resonance vector
is positive, so scaling \(x\) can set the coefficient to one.
The residual equation \(\det JF=1\) is the coefficient ideal of
this chart over \(\mathbb Q\).

Three weights use smaller presentations of their charts.

- For \((1,-1,-2)\), use the
  [four ordinary sextic charts already proved](SEXTIC_KELLER_GRADED_112.md).
- For \((1,-1,-1)\), put \(u=xy,\ v=xz\). The normalized map is
  \((x\Lambda,B/x,C/x)\), where \(\Lambda(0)=1,\ \deg\Lambda\le2\),
  \(B=u+O((u,v)^2),\ C=v+O((u,v)^2)\), and
  \(\deg B,\deg C\le3\) in \(u,v\). The quotient coordinates are
  \(\Lambda B,\Lambda C\), and the determinant condition is
  \[
   \Lambda\{B,C\}+B\{\Lambda,C\}+C\{B,\Lambda\}=1.
  \]
  This follows by the chain rule in \((x,u,v)\), whose source
  coordinate change has Jacobian \(x^2\). A linear change in
  \(y,z\) gives the same change of \(u,v\), and inverse target
  normalization preserves the identity first derivative. If the linear
  part of \(\Lambda\) is nonzero, take it as \(u\).
  Completing the square in \(v\) if its square coefficient is
  nonzero, and otherwise using a shear and scaling if the \(uv\)
  coefficient is nonzero, gives \(1+u+v^2+c u^2\),
  \(1+u+uv\), or \(1+u+c u^2\). If the linear part is zero,
  diagonalizing the nonzero binary quadratic gives \(1+uv\) or
  \(1+u^2\). These are the five balanced charts.
- For \((1,0,-1)\), consider the positive and negative multipliers
  \(F_1/x\) and \(F_3/z\), as polynomials in \(xz,y\).
  If either coefficient of \(xz\) is nonconstant in \(y\),
  exchange \(x,z\) and reverse the parameter if needed, so this
  is true in the positive multiplier. Translate the fixed-axis
  variable \(y\) to a root of that coefficient and renormalize
  the target by the corresponding translation. The first derivative
  remains the identity by the fixed-axis calculation. The term
  \(x^2z\) in the positive component now has coefficient zero.
  Thus the first chart, where that term has normalized coefficient
  one, is only needed when the coefficients of \(xz\) in both
  multipliers are constant in \(y\). The other charts use the
  full support from the fixed-axis calculation.

The [exact checker](verify_sextic_graded_action.py) generates the
resonance pairs, chart supports and determinant coefficients. It
replays 349 charts with no zero weights besides the two special
weights, 20 mixed-sign charts with a fixed axis, and the five balanced
charts. Each ideal has reduced Gröbner basis \(\{1\}\) over
\(\mathbb Q\) in Singular, hence the chart has no point over any
characteristic-zero field. The four omitted charts are the ones
for the known example's weights. These calculations exclude all
nonlinear positive-weight outputs in the mixed-sign case.

The replay command is `.venv/bin/python research/jacobian/monodromy-theory/verify_sextic_graded_action.py`
from the repository root, using SymPy 1.14.0 and Singular 4.3.2.
The [checker for the known example's weights](verify_sextic_graded_112.py)
is run separately. Both results were checked internally.

Inversion descends to the original field by uniqueness of the
normalized formal inverse. The rational degree-seven example from
[Gao, Theorem 3.3](sources/gao-2608.00222v2.pdf) gives the
upper bound for every characteristic-zero field.
