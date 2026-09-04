# Smooth-quartic squarefree-line generic gate

## Status

This note proves `HC4NHM16`.  It treats one of the four basepoint-free
simple-residual-line rows isolated by `HC4NHM14`, namely

\[
d_0=(x^2,y^2,0).
\]

The generic point of this incidence is empty: the exact reciprocal and
Hessian equations force the quadratic matrix to have zero determinant.
This is a generic-stratum theorem, not a closure of the complete row.  A
proper Zariski-closed exceptional parameter locus, the complementary
residual-line chart, and the other three basepoint-free rows remain.

Replay the exact characteristic-zero calculation with

```bash
.venv/bin/python scripts/verify_hc4_smooth_quartic_squarefree_line_generic.py
```

The checker uses SymPy to construct the coefficient equations and Singular
4.4.1 for the rational-function-field standard basis.

For repository maintenance, `--audit-existing-only` verifies the exact hash
of the imported 81-equation builder without constructing the system or
invoking Singular.  It is a provenance check, not a replacement proof.

## 1. The squarefree-line chart

Normalize the reciprocal line to \(z=0\).  In the squarefree binary-cubic
row with zero normal class, put

\[
s_3=\frac{x^3+y^3}{3}+z^2(ux+vy)+wz^3,
\qquad d=\nabla s_3,
\tag{1.1}
\]

and use the complete boundary matrix from (4.2) of `HC4NHM14`, scaled so
that \(\lambda=-1\):

\[
A_0=\begin{pmatrix}
0&0&-y^2\\
0&0&x^2\\
-y^2&x^2&p x^2+qxy+r y^2
\end{pmatrix}.
\tag{1.2}
\]

Write \(A=A_0+zB\), where \(B\) is a general symmetric matrix of linear
forms.  Its six upper-triangular entries, in the order

\[
(00),(01),(02),(11),(12),(22),
\]

have coefficient triples

\[
(b_0,b_1,b_2),\ldots,(b_{15},b_{16},b_{17}).
\tag{1.3}
\]

On the affine residual-line chart with nonzero \(y\)-coefficient, write

\[
\ell=y+\tau x+\sigma z.
\tag{1.4}
\]

Define the exact polynomial quotients

\[
C=\frac{\operatorname{adj}(A)+dd^{\mathsf T}}{z},\qquad
e=\frac{Ad}{z},\qquad
R=\frac{\det A}{z}.
\tag{1.5}
\]

Their divisibility by \(z\) follows from the boundary identities in
`HC4NHM14` and is checked during construction.

## 2. Necessary reciprocal equations

Every smooth-quartic solution in this chart obeys three groups of equations:

1. the nine Hessian-curl identities for the symmetric cubic matrix \(C\);
2. \(R|_{y=-\tau x-\sigma z}=0\), equivalently \(\ell\mid R\);
3. the scalar compatibility

   \[
   \ell\bigl(R+e^{\mathsf T}d\bigr)=zR.
   \tag{2.1}
   \]

Indeed, if \(R=Q\ell\), then (2.1) is exactly
\(d^{\mathsf T}e=Qa\) with \(z=\ell+a\), because \(\lambda=-1\).
Thus this system is a relaxation of the original reciprocal problem, and
showing that it has only \(\det A=0\) solutions is enough for exclusion.

Coefficient extraction gives 81 distinct nonzero equations: 10 have total
deformation degree one, 52 have degree two, and 19 have degree three.

## 3. Generic standard basis

Work over

\[
K=\mathbb Q(\tau,\sigma,p,q,r,b_{15},b_{16},b_{17})
\tag{3.1}
\]

and take

\[
b_0,\ldots,b_{14},u,v,w
\tag{3.2}
\]

as polynomial variables.  The ten purely linear equations have rank ten
over \(K\).  A visible common pivot is

\[
\begin{aligned}
\Delta={}&(3p^2-qr)\tau^5+(9pr-q^2)\tau^4
 +(18r^2-6pq)\tau^3\\
&+(18p^2-6qr)\tau^2+(9pr-q^2)\tau+(3r^2-pq),
\end{aligned}
\tag{3.3}
\]

which is nonzero as a polynomial.  Reducing the remaining 71 equations by
that linear basis and applying `slimgb` gives

\[
(w,v,u,b_{14},b_{11},b_8,b_5,b_2).
\tag{3.4}
\]

Combining (3.4) with the ten linear pivots gives the reduced basis

\[
(w,v,u,b_{14},b_{13},\ldots,b_0).
\tag{3.5}
\]

Thus the generic relaxed fiber has

\[
b_0=\cdots=b_{14}=u=v=w=0.
\tag{3.6}
\]

Only the bottom-right normal coefficients remain, so

\[
A=\begin{pmatrix}
0&0&-y^2\\
0&0&x^2\\
-y^2&x^2&p x^2+qxy+r y^2+z(b_{15}x+b_{16}y+b_{17}z)
\end{pmatrix},
\tag{3.7}
\]

and \(\det A=0\).  Hence \(R=0\) and the putative quartic quotient \(Q\)
is zero, contradicting the smooth nonzero quartic target.

## 4. Result and exact remainder

> **Theorem `HC4NHM16` -- Generic squarefree-line gate.**  In the
> basepoint-free simple-residual-line type
> \(d_0=(x^2,y^2,0)\) of `HC4NHM14`, the generic point of the affine
> incidence chart \(\ell=y+\tau x+\sigma z\) admits no reciprocal Hessian
> solution with nonzero quartic quotient.  Over the total parameter
> function field the 81 necessary coefficient equations force (3.6), after
> which \(\det A=0\).

The word *generic* is essential.  A function-field standard basis proves
emptiness over some nonempty Zariski-open subset; it does not classify the
specializations where coefficients inverted during the calculation vanish.
Equation (3.3) records the first visible exceptional divisor, but is not
claimed to be the complete denominator of the certificate.

The subsequent invariant identification `HC4NHM20` proves that (3.3) is the
first polar of the binary resultant `Res(s^3+t^3,H)` in an explicit
quadratic direction `K_tau`.  Its generic `(p,q,r)`-fiber is a smooth conic;
at the fifteen squarefree roots of `Res(s^3+t^3,K_tau)` it splits into one
resultant line and one residual polar line.  This is not itself the
tangent/bitangent/flex discriminant of the prospective quartic `Q`; see
[`HC4_SMOOTH_QUARTIC_PIVOT_POLAR_GEOMETRY.md`](HC4_SMOOTH_QUARTIC_PIVOT_POLAR_GEOMETRY.md).

`HC4NHM17` enters nine exact pieces of that first visible divisor.  It
closes the central quadratic-zero locus generically and at its first
algebraic pivot, the generic charts and first algebraic pivot of the fiber
\(\tau=0\), and both
components of the fiber \(\tau=-1\) together with their first displayed
secondary pivots.  See
[`HC4_SMOOTH_QUARTIC_SQUAREFREE_LINE_EXCEPTIONAL_SLICES.md`](HC4_SMOOTH_QUARTIC_SQUAREFREE_LINE_EXCEPTIONAL_SLICES.md).

The generic point of the complete polar conic is subsequently excluded by
`HC4NHM22`.  Its universal point `[1:3:1]` gives a rational parametrization,
and an eight-coefficient radical certificate puts the sixth power of every
active deformation coordinate in the necessary ideal.  See
[`HC4_SMOOTH_QUARTIC_POLAR_CONIC_GENERIC_GATE.md`](HC4_SMOOTH_QUARTIC_POLAR_CONIC_GENERIC_GATE.md).

The exact remaining work inside action 1 is therefore:

1. continue the lower denominator strata left by the generic polar-conic
   exclusion `HC4NHM22`, the non-generic two-line fibers, and the
   complementary residual-line chart;
2. perform the corresponding generic and exceptional eliminations for
   \((x^2,y^2,xy)\), \((2xy,x^2,y^2)\), and \((x^2,0,y^2)\);
3. only afterward move to the five basepointed simple-line rows with their
   forced contact allocations.

The doubled-line rows and the scalar-degenerate \(\mu=0\) row are separate
actions and are not treated here.
