# The foundational incidence construction

This note contains only the projective geometry of the foundational
factorization construction: arbitrary cubic hyperplanes, the three contact
orbits with the twisted cubic, the Grothendieck-class obstruction, and the
exceptional tangent-nonosculating orbit.  The normalized complete
intersection, its coordinates and inverse, and its étaleness are proved once
in the [normalized factorization model](NORMALIZED_FACTORIZATION_MODEL.md).

Work over a characteristic-zero field `k`.  Put

\[
 V_i=\operatorname{Sym}^i(k^2)
\]

and let `R(L,Q)=Res(L,Q)` for `(L,Q) in V_1 x V_2`.  For a nonzero
`ell in V_3^*`, define

\[
 X_\ell=\{(L,Q):R(L,Q)=1,\ \ell(LQ)=1\}
\]

and

\[
 U_\ell=
 \bigl(\mathbb P(V_1)\times\mathbb P(V_2)\bigr)
 \setminus\bigl(\{R=0\}\cup\{\ell(LQ)=0\}\bigr).       \tag{1}
\]

The algebraic model proves that multiplication
`X_ell -> {C:ell(C)=1}` is étale for every `ell`.

## Projective normalization and generic degree

Projectivization gives a canonical isomorphism

\[
 \boxed{X_\ell\xrightarrow{\sim}U_\ell.}              \tag{2}
\]

Indeed, choose representatives of a point of `U_ell` and put

\[
 m=\ell(LQ),\qquad r=R(L,Q).
\]

Under `(L,Q)->(alpha L,beta Q)`, these values become

\[
 m\longmapsto\alpha\beta m,\qquad
 r\longmapsto\alpha^2\beta r.
\]

There is a unique rescaling for which both normalized values equal one:

\[
 \boxed{\alpha={m\over r},\qquad \beta={r\over m^2}.} \tag{3}
\]

If the initial representatives are changed, the scalars in (3) change by
the inverse representative scalars, so the normalized pair is unchanged.
The formulas are regular on `U_ell` and invert projectivization.

For a squarefree cubic `C`, a point above `C` is the choice of one of its
three linear factors; (3) then fixes the two factor scalings uniquely.
Consequently normalized multiplication is generically three-to-one.  This
projective marking explains the degree, but not why one particular source is
affine space.

## The three hyperplane orbits

Assume now that `k` is algebraically closed.  The restriction of `ell` to the
twisted cubic is a binary cubic.  Its zero divisor has exactly one of the
three multiplicity partitions

\[
 (1,1,1),\qquad(2,1),\qquad(3).                      \tag{4}
\]

These are precisely the projective hyperplane orbits.  Indeed, the zero
divisor is invariant under projective change of variables, while `PGL_2` is
transitive on ordered triples of distinct points and on ordered pairs of
distinct points.  Projectively the `SL_2` and `PGL_2` orbits agree.

The corresponding normalized sources have classes

\[
\begin{array}{c|c}
\text{contact type}&[X_\ell]\\ \hline
(1,1,1)&\mathbb L^3-\mathbb L\\
(2,1)&\mathbb L^3\\
(3)&\mathbb L^3-\mathbb L^2.
\end{array}                                          \tag{5}
\]

Over `C`, the Hodge--Deligne realization therefore gives

\[
 \boxed{X_\ell\simeq\mathbb A^3
 \quad\Longleftrightarrow\quad
 \ell\text{ has contact type }(2,1).}               \tag{6}
\]

The positive orbit is tangent but nonosculating.

## Grothendieck-class calculation

Project (1) to the marked linear factor
`[L] in P(V_1)=P^1`.  In the residual `P^2`, the resultant condition removes
the line of quadratics divisible by `L`, and the hyperplane condition removes
the kernel of

\[
 Q\longmapsto\ell(LQ).                               \tag{7}
\]

Away from a multiple contact point the two removed lines are distinct, so
the fiber has class

\[
 [\mathbb P^2]-2[\mathbb P^1]+[\mathrm{pt}]
 =\mathbb L^2-\mathbb L.                             \tag{8}
\]

At a double contact point, (7) vanishes precisely on the resultant line.  The
two removed lines coincide and the fiber is `A^2`, of class `L^2`.  At a
triple contact point, (7) is identically zero, so the pulled-back hyperplane
contains the whole residual fiber and the open fiber is empty.

For type `(1,1,1)`, every marked point has the first fiber type, hence

\[
 [X_\ell]=(1+\mathbb L)(\mathbb L^2-\mathbb L)
 =\mathbb L^3-\mathbb L.                             \tag{9}
\]

For type `(2,1)`, the double point contributes `L^2`, while its affine-line
complement in the base contributes the ordinary fibers:

\[
 [X_\ell]=\mathbb L^2+
 \mathbb L(\mathbb L^2-\mathbb L)=\mathbb L^3.       \tag{10}
\]

For type `(3)`, the triple-point fiber is empty and the remaining base is
`A^1`, so

\[
 [X_\ell]=\mathbb L(\mathbb L^2-\mathbb L)
 =\mathbb L^3-\mathbb L^2.                           \tag{11}
\]

This proves (5) and obstructs affine three-space in the first and third
orbits over `C`.

## Why type `(2,1)` is exceptional

Suppose the contact divisor is `2p+q` with `p!=q`.  The group `PGL_2` is
transitive on such ordered pairs, and binary-form multiplication commutes
with its action.  A projective change of variables may multiply the resultant
and hyperplane functional by nonzero scalars `delta` and `eta`; the additional
factor rescaling

\[
 L\longmapsto {\eta\over\delta}L,
 \qquad Q\longmapsto {\delta\over\eta^2}Q
\]

restores both normalized equations.  Hence every `(2,1)` source is isomorphic
to the representative

\[
 a^2e-abd+b^2c=1,\qquad ad+bc=1.                    \tag{12}
\]

Proposition 1 of the
[normalized factorization model](NORMALIZED_FACTORIZATION_MODEL.md) gives a
global polynomial isomorphism from (12) to `A^3`, including the divisor
`a=0`.  That same note owns the explicit forward coordinates, polynomial
inverse, étaleness, determinant, and linear equivalence with the foundational
map.

Thus the construction separates into three independent mechanisms:

\[
\begin{array}{rcl}
\text{étaleness}&\Longleftarrow&
 \text{coefficient--resultant algebra for every }\ell,\\
\text{generic degree three}&\Longleftarrow&
 \text{three choices of a marked linear factor},\\
X_\ell\simeq\mathbb A^3&\Longleftarrow&
 \text{the exceptional }(2,1)\text{ hyperplane orbit}.
\end{array}
\]

The orbit and class assertions are internal repository results with no
recorded external review.  Direct finite-field checks of the three class
formulas are included in
[`verify_normalized_factorization_slice.py`](../scripts/verify_normalized_factorization_slice.py).
