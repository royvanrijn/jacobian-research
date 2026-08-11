# Smooth-quartic reciprocal-factorization frontend

## Status

This note proves `HC4NHM14`.  It continues the clean generic-corank-one
quartic-denominator branch of `HC4NHM1` when the minimal denominator is a
smooth irreducible quartic.  It converts the cubic Hessian factorization into
a reciprocal quadratic factorization and reduces its residual-line boundary
to a finite list of gradient-kernel types.  It is a frontend theorem, not an
exclusion of the smooth-quartic packet.

Replay the exact matrix identities and boundary representatives with

```bash
.venv/bin/python scripts/verify_hc4_smooth_quartic_reciprocal_frontend.py
```

The checker verifies the polynomial matrix identities, the displayed binary
representatives, their common-factor degrees, and the four basepoint-free
quadratic-matrix families.  Exhaustiveness of the binary-cubic root split is
the elementary written argument in Section 3.

## 1. The clean smooth-quartic data

Let

\[
 C=\operatorname{Hess}(h_5),\qquad
 \det C=Q^2\ell,
\tag{1.1}
\]

where \(Q\) is a smooth irreducible ternary quartic and \(\ell\) is a
nonzero linear form.  The clean module section of `HC4NHM1` consists of a
cubic vector \(e\), a quadratic gradient \(d=\nabla s_3\), and a linear form
\(a\) satisfying

\[
 Ce=Qd,\qquad d^{\mathsf T}e=Qa.
\tag{1.2}
\]

On \(Q=0\), the vector \(e\) is a nowhere-vanishing generator of the
kernel line.  The two symmetric rank-one adjugates on the smooth projective
curve differ by a nonzero scalar, so for a symmetric quadratic matrix \(A\)
and \(\lambda\in K^\times\),

\[
 \operatorname{adj}(C)=\lambda ee^{\mathsf T}+QA.
\tag{1.3}
\]

Multiplying (1.3) by \(C\) on either side gives the paired equations

\[
 \boxed{
 \begin{aligned}
 CA+\lambda de^{\mathsf T}&=Q\ell I,\\
 AC+\lambda ed^{\mathsf T}&=Q\ell I.
 \end{aligned}}
\tag{1.4}
\]

Put

\[
 \mu=\ell-\lambda a.
\tag{1.5}
\]

Multiplication of the first equation in (1.4) by \(d\), followed by (1.2),
gives

\[
 C(Ad-\mu e)=0.
\]

Since \(\det C\ne0\),

\[
 \boxed{Ad=\mu e.}
\tag{1.6}
\]

## 2. The reciprocal quadratic factorization

Assume first that \(\mu\ne0\).  The rank-one determinant identity applied
to (1.3) is

\[
 Q\det A+\lambda e^{\mathsf T}\operatorname{adj}(A)e
 =Q^2\ell^2.
\tag{2.1}
\]

Equation (1.6) implies

\[
 \mu e^{\mathsf T}\operatorname{adj}(A)e
 =\det(A)e^{\mathsf T}d
 =Q a\det A.
\tag{2.2}
\]

Multiply (2.1) by \(\mu\), use (2.2), and use
\(\mu+\lambda a=\ell\).  Cancellation in the polynomial domain gives

\[
 \boxed{\det A=Q\ell\mu.}
\tag{2.3}
\]

Moreover,

\[
 A(\mu C+\lambda dd^{\mathsf T})
 =\mu(AC+\lambda ed^{\mathsf T})
 =Q\ell\mu I.
\]

Together with (2.3), this proves

\[
 \boxed{\operatorname{adj}(A)=\mu C+\lambda dd^{\mathsf T}.}
\tag{2.4}
\]

In particular,

\[
 \operatorname{adj}(A)e=Q\ell d,
 \qquad
 C=\frac{\operatorname{adj}(A)-\lambda dd^{\mathsf T}}{\mu}.
\tag{2.5}
\]

Thus the original \(3\)-by-\(3\) cubic Hessian problem has acquired a
reciprocal \(3\)-by-\(3\) quadratic determinantal problem.  The gradient
condition on \(d\) and the Hessian integrability of the quotient in (2.5)
remain essential.

If \(\mu=0\), then (1.6) gives \(Ad=0\), and \(d\ne0\) over the fraction
field by (1.2); hence \(\det A=0\).  This scalar-degenerate branch does not
carry the residual-line factorization (2.3) and must be treated separately.

## 3. The residual-line boundary

Suppose \(\mu\ne0\), normalize \(\mu=z\), and put \(L=(z=0)\).  Write

\[
 s_3=F_3(x,y)+zG_2(x,y)+z^2H_1(x,y)+cz^3,
\tag{3.1}
\]

so that

\[
 d_0=d|_L=(F_x,F_y,G).
\tag{3.2}
\]

Writing \(A=A_0+zA_1\), restriction of (1.6) and (2.4) gives

\[
 A_0d_0=0,
 \qquad
 \operatorname{adj}(A_0)=\lambda d_0d_0^{\mathsf T}.
\tag{3.3}
\]

The first normal determinant coefficient is therefore

\[
 \left.\frac{\det A}{z}\right|_L
 =\operatorname{tr}(\operatorname{adj}(A_0)A_1)
 =\lambda d_0^{\mathsf T}A_1d_0.
\tag{3.4}
\]

If \(\ell\) and \(L\) are distinct, (2.3) and (3.4) give

\[
 \boxed{
 (Q|_L)(\ell|_L)=\lambda d_0^{\mathsf T}A_1d_0.
 }
\tag{3.5}
\]

Consequently, for \(g=\gcd(F_x,F_y,G)\),

\[
 \boxed{g^2\mid (Q|_L)(\ell|_L).}
\tag{3.6}
\]

If instead \(\ell\) is proportional to \(z\), the residual line is doubled
in \(\det A\), and the first normal condition is

\[
 \boxed{d_0^{\mathsf T}A_1d_0=0.}
\tag{3.7}
\]

### 3.1 Complete binary-gradient list

Changes of \(x,y\), scaling of \(z\), and the shears
\(x\mapsto x+\alpha z,\ y\mapsto y+\beta z\) act on (3.1) by the usual
binary change on \(F\) and by

\[
 G\longmapsto G+\alpha F_x+\beta F_y.
\tag{3.8}
\]

The root type of a nonzero binary cubic is squarefree, double, or triple.
Splitting \(G\) modulo the span of its two derivatives gives the following
complete list.  The final zero row occurs only on the doubled residual-line
boundary: on a simple residual line it would make \(A_0\) generically rank
at most one, while \(\det A\) has generic order one.

| binary-cubic type | normal class | \(d_0=(F_x,F_y,G)\) | \(\deg g\) |
|---|---|---|---:|
| squarefree | outside derivative span | \((x^2,y^2,xy)\) | 0 |
| squarefree | in derivative span | \((x^2,y^2,0)\) | 0 |
| double root | outside derivative span | \((2xy,x^2,y^2)\) | 0 |
| double root | in derivative span | \((2xy,x^2,0)\) | 1 |
| triple root | square normal class | \((x^2,0,y^2)\) | 0 |
| triple root | simple normal class | \((x^2,0,xy)\) | 1 |
| triple root | zero normal class | \((x^2,0,0)\) | 2 |
| \(F=0\) | squarefree \(G\) | \((0,0,xy)\) | 2 |
| \(F=0\) | double-root \(G\) | \((0,0,x^2)\) | 2 |
| \(F=G=0\) | zero | \((0,0,0)\) | -- |

For example, a degree-one common factor in the distinct-line case forces
the residual line to meet \(Q\cup(\ell=0)\) with total multiplicity at
least two at that point.  A squarefree degree-two common factor forces a
bitangent-type contact allocation, while a double degree-two common factor
forces flex or hyperflex contact after accounting for \(\ell|_L\).

## 4. Basepoint-free quadratic boundary matrices

The four rows with \(\deg g=0\) admit small complete \(A_0\)-families.
They are the natural first exact eliminations because no contact divisor is
available to absorb a kernel basepoint.

For \(d_0=(x^2,y^2,xy)\),

\[
A_0=\begin{pmatrix}
p y^2&-qxy&y(qy-px)\\
-qxy&r x^2&x(qx-ry)\\
y(qy-px)&x(qx-ry)&p x^2-2qxy+r y^2
\end{pmatrix},
\quad \lambda=pr-q^2.
\tag{4.1}
\]

For \(d_0=(x^2,y^2,0)\),

\[
A_0=\begin{pmatrix}
0&0&-t y^2\\
0&0&t x^2\\
-t y^2&t x^2&H_2
\end{pmatrix},
\quad \lambda=-t^2,
\tag{4.2}
\]

where \(H_2\) is an arbitrary binary quadratic.  For
\(d_0=(2xy,x^2,y^2)\),

\[
A_0=\begin{pmatrix}
(p x^2+2qxy+r y^2)/4&-y(px+qy)/2&-x(qx+ry)/2\\
-y(px+qy)/2&p y^2&qxy\\
-x(qx+ry)/2&qxy&r x^2
\end{pmatrix},
\quad \lambda=(pr-q^2)/4.
\tag{4.3}
\]

Finally, for \(d_0=(x^2,0,y^2)\),

\[
A_0=\begin{pmatrix}
0&-t y^2&0\\
-t y^2&H_2&t x^2\\
0&t x^2&0
\end{pmatrix},
\quad \lambda=-t^2.
\tag{4.4}
\]

Direct multiplication gives \(A_0d_0=0\) and
\(\operatorname{adj}(A_0)=\lambda d_0d_0^{\mathsf T}\) in all four
families.  Generic boundary rank requires \(pr-q^2\ne0\) in (4.1),
(4.3), and \(t\ne0\) in (4.2), (4.4).

## 5. Result and next rows

> **Theorem `HC4NHM14` -- Smooth-quartic reciprocal frontend.**  In the
> clean generic-corank-one smooth-quartic denominator packet of `HC4NHM1`,
> put \(\mu=\ell-\lambda a\).  If \(\mu=0\), the reciprocal quadratic
> matrix is generically singular and satisfies \(Ad=0\).  If
> \(\mu\ne0\), then
> 
> \[
> \det A=Q\ell\mu,qquad
> \operatorname{adj}(A)=\mu C+\lambda dd^{\mathsf T},qquad
> Ad=\mu e.
> \]
> 
> On the residual line \(\mu=0\), the gradient kernel has exactly the ten
> types in Section 3.1.  All ten are a priori doubled-line types, with the
> zero type as the extra degeneration.  On the distinct-line boundary the
> zero type is impossible and every common factor \(g\) obeys
> \(g^2\mid(Q|_L)(\ell|_L)\).  Hence the smooth-quartic packet splits into
> one scalar-degenerate row, nine simple residual-line rows, and ten doubled
> residual-line rows.

The next calculation begins with the four basepoint-free simple-line rows
(4.1)--(4.4).  Extend \(A_0\) by a general \(zA_1\), form

\[
C=\frac{\operatorname{adj}(A)-\lambda dd^{\mathsf T}}{z},
\]

impose Hessian curl equations, and require
\(\det A/(z\ell)\) to be a smooth irreducible quartic.  `HC4NHM16` carries
this out at the generic point of the squarefree-line row (4.2) and excludes
that generic stratum.  `HC4NHM17` further excludes nine exact slices of its
first visible exceptional divisor.  The rest of that divisor and the
complementary residual-line chart remain, as do the other three
basepoint-free families.
The five basepointed simple-line rows should then be treated on their forced
tangent/bitangent/flex incidence strata.  The doubled-line and \(\mu=0\)
branches remain separate.  See
[`HC4_SMOOTH_QUARTIC_SQUAREFREE_LINE_GENERIC_GATE.md`](HC4_SMOOTH_QUARTIC_SQUAREFREE_LINE_GENERIC_GATE.md).
