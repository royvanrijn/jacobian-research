# Controlled-boundary suspensions

This note formalizes the common determinant mechanism behind the weighted
and cancellation constructions, proves the elementary plane-core normal
forms that can already be classified, and isolates the first obstruction to
a genuinely independent two-boundary cancellation chart.  It is a scoped
starting point, not a classification of all polynomial Keller maps.

Work over a characteristic-zero field `k`.

## 1. A suspension is a square, not only a composition

The weighted and cancellation constructions fit a common commutative square

```text
 X = A^3  -------- F -------->  Y = A^3
    |                              |
    | alpha                        | beta
    v                              v
    Z  -------- Phi -------->      T
```

where `alpha` and `beta` are dominant rational maps between smooth
threefolds, `Phi` is the product of a plane core with a parameter coordinate,
and

\[
 \beta\circ F=\Phi\circ\alpha.                         \tag{1}
\]

Let `E=(D=0)` be the sole nonunit divisor in the core Jacobian and suppose

\[
 \det D\Phi=uD^r,
 \qquad u\in k^*.                                     \tag{2}
\]

Taking determinants in (1) gives

\[
 (\det D\beta\circ F)\det DF
 =u(D\circ\alpha)^r\det D\alpha.                       \tag{3}
\]

If `F` has constant nonzero Jacobian, the following equality of principal
Weil divisors holds:

\[
 \boxed{
 \operatorname{div}(\det D\alpha)
 +r\operatorname{div}(D\circ\alpha)
 =F^*\operatorname{div}(\det D\beta).}                 \tag{4}
\]

We call (4) the **determinant ledger**.  Pullback through a rational chart is
understood as the divisor of the pulled-back rational function, so its
coefficients may be negative.

Conversely, equality of the ledger says that `det DF` is a global unit.  It
says that this unit is constant when `O(X)^*=k^*`, in particular for
`X=A^3`.  On an arbitrary smooth threefold with nonconstant units, (4) alone
does not strengthen “global unit” to “constant Jacobian.”

For cancellation, `beta=id`, `D\circ\alpha=A^{-1}`, and
`det D alpha=-A^r`; the two terms on the left of (4) cancel.  For the
weighted square,

\[
 \det D\alpha=b_0x^3\gamma^2,
 \quad D\circ\alpha=\gamma,
 \quad \det D\beta=-cC^3,
 \quad C=x\gamma,
\]

so (4) is the identity

\[
 3\operatorname{div}(x)+3\operatorname{div}(\gamma)
 =3\operatorname{div}(C).
\]

Consequently a datum `(Phi,D,r,alpha)` describes cancellation, but a common
definition covering the weighted construction must also retain the target
chart `beta`.

## 2. Coordinate-preserving plane cores

Let

\[
 \phi:\mathbb A^2_{w,q}\longrightarrow\mathbb A^2_{q,t},
 \qquad \phi(w,q)=(q,T(w,q)).                           \tag{5}
\]

Then

\[
 \det D\phi=-T_w.                                      \tag{6}
\]

This gives the following elementary but useful normal-form lemma.

### Proposition 2.1 -- simple section normal form

Suppose the critical divisor is the reduced section

\[
 E=(q-h(w)=0)
\]

and

\[
 \det D\phi=u(h(w)-q),\qquad u\in k^*.
\]

Choose `H` with `H'=h`.  Up to a target shear preserving `q` and a nonzero
scaling of `t`, the map is

\[
 \boxed{(w,q)\longmapsto(q,wq-H(w)).}                  \tag{7}
\]

Indeed (6) gives `T_w=u(q-h(w))`, hence

\[
 T=u(wq-H(w))+g(q).
\]

The target automorphism
`(q,t)\mapsto(q,u^{-1}(t-g(q)))` gives (7).  Thus the
universal marked-root incidence is forced by four hypotheses: a preserved
coordinate, a critical divisor which is a section, reduced ramification,
and a constant Jacobian unit along that section.

### Proposition 2.2 -- primitive form for an arbitrary boundary power

If instead

\[
 \det D\phi=uD(w,q)^r,
\]

then, up to a target shear preserving `q`,

\[
 T(w,q)=-u\int_0^wD(v,q)^r\,dv.                        \tag{8}
\]

Formula (8) classifies coordinate-preserving cores once `D` is fixed, but
does not classify the divisors `D` that admit a polynomial Keller suspension.
That suspension condition is the restrictive part of the problem.

## 3. The first intrinsic split

For the weighted core, `D=q-H'(w)`, so `E` is a graph and

\[
 E\simeq\mathbb A^1,
 \qquad \mathcal O(E)^*=k^*.
\]

For the cancellation core at fixed `P`, put

\[
 D=1-s(Q-Ps)^m,
 \qquad Y=Q-Ps.
\]

On `D=0` one has `sY^m=1`, and conversely

\[
 s=Y^{-m},\qquad Q=Y+PY^{-m}.
\]

Therefore

\[
 k[E]=k[Y,Y^{-1}],
 \qquad E\simeq\mathbb G_m.                            \tag{9}
\]

The affine type, unit rank, and number of places at infinity of the critical
normalization separate the two plane cores before any threefold boundary
calculation.  The ramification index at the generic discriminant point is
`r+1`; the two pole vectors in (9) retain the cancellation exponent `m`.

## 4. A first independent two-boundary ansatz

The current cancellation chart has one source valuation
`A=1+xf(y)`.  To test the smallest genuine enlargement, take two distinct
functions `f_1,f_2` for which

\[
 A_1=1+xf_1(y),\qquad A_2=1+xf_2(y)
\]

are algebraically independent.  Consider the multiplicative triangular
chart

\[
 B=A_1^{b_1}A_2^{b_2}z+g(x,y),
 \quad P=A_1^{a_1}A_2^{a_2}B,
\]

\[
 s=xA_1^{d_1}A_2^{d_2},
 \qquad Q=y+sP.                                       \tag{10}
\]

This is the direct two-divisor analogue of the one-inverse-variable
triangular reconstruction skeleton.

### Proposition 4.1 -- third-divisor obstruction

The chart Jacobian is

\[
 \det\frac{\partial(s,P,Q)}{\partial(x,y,z)}
 =-A_1^{a_1+b_1+d_1-1}A_2^{a_2+b_2+d_2-1}N_{d_1,d_2}, \tag{11}
\]

where

\[
 N_{d_1,d_2}
 =(1+d_1+d_2)A_1A_2-d_1A_2-d_2A_1.                   \tag{12}
\]

The polynomial `N_(d_1,d_2)` is a monomial in `A_1,A_2` exactly for

\[
 (d_1,d_2)=(0,0),\quad(-1,0),\quad(0,-1).             \tag{13}
\]

To prove (11), replace the final coordinate `Q=y+sP` by `y`, which is a
target-triangular determinant-one operation.  Then the determinant is
`-P_z s_x`, and differentiating `s` gives (11)-(12).  Since (12) is supported
on the three distinct monomials `A_1A_2,A_1,A_2`, it is a monomial only when
two of their coefficients vanish, giving exactly (13).

The first case in (13) has no inverse factor in `s`; each of the other two
inverts only one of `A_1,A_2`.  If both exponents are genuinely active, the
Jacobian acquires the additional divisor `N_(d_1,d_2)=0`.  Hence this
smallest shared-reconstruction-variable ansatz cannot produce a clean
two-boundary determinant ledger.  A new construction must do at least one of
the following:

1. accept and control a third boundary divisor;
2. use two reconstruction variables;
3. use a nonmultiplicative source chart;
4. use a nontrivial target ledger to absorb `N_(d_1,d_2)`.

This proposition is an exact obstruction for (10), not evidence that all
two-boundary suspensions are impossible.

## 5. A realistic classification target

The first attainable theorem is an **elementary one-boundary suspension
dichotomy** under the following explicit restrictions:

1. the plane core preserves a coordinate;
2. the critical divisor is smooth, rational, and has at most two places at
   infinity;
3. the suspension square is divisor-minimal;
4. reconstruction uses one primitive variable;
5. the vertical charts are either polynomial weighted charts or rank-one
   triangular affine modifications.

Under these hypotheses, Proposition 2.1 supplies the weighted core when the
critical curve is an affine line, while the rigidity theorem for the current
cancellation ansatz supplies the punctured-line branch.  The missing global
step is a classification of divisor-minimal one-valuation rational charts of
`A^3` as rank-one affine modifications.  Without that step, exhaustiveness
must not be claimed.

## 6. Where the present method stops

The next concrete search should allow a target chart in (10) or add a second
primitive reconstruction variable and solve the resulting mixed finite-jet
conditions.  In contrast, immediately attempting to classify arbitrary
birational charts of `A^3` with prescribed Jacobian divisor is an open-ended
affine-modification problem.  That is the point at which this pathway reaches
diminishing returns unless a new structural theorem or a surviving low-degree
example supplies additional leverage.

The first target-ledger example now exists after leaving affine space:
the [three-linear-factor Cox ledger](../extended-geometry/COX_LEDGER_THREE_FACTOR.md)
keeps its rank-one boundary unit and cancels it with one primitive coordinate,
giving a constant-Jacobian finite etale morphism between affine fourfold
boundary complements.  Its reciprocal coordinate does not extend across the
collision divisor, so it does not solve the polynomial `A^4` chart problem.
The [all-arity extension](../extended-geometry/COX_LEDGER_LINEAR_FACTORS.md)
constructs the corresponding separated ledger for every `s>=3` and shows
that the unit-rank dimension bound is false without the separate-character
reconstruction requirement.
In the first arity, the
[oriented cubic target chart](../extended-geometry/ORIENTED_CUBIC_COX_CHART.md)
absorbs the primitive discriminant character polynomially, extends across
one collision divisor, and leaves the other two ordered collision branches
as distinct dicritical divisors.

The exact determinant identities and bounded exponent audit are checked by
[`verify_controlled_boundary_suspensions.py`](../scripts/verify_controlled_boundary_suspensions.py).
The later
[minimal-boundary classification note](MINIMAL_BOUNDARY_CLASSIFICATION.md)
incorporates the saturated-link, boundary-monotonicity, and reciprocal-branch
advances and is now the canonical statement of the open classification
problem.
