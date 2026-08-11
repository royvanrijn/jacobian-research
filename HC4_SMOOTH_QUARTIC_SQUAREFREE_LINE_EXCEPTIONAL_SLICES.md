# Smooth-quartic squarefree-line exceptional slices

## Status

This note proves `HC4NHM17`.  It continues the first, function-field
squarefree-line exclusion `HC4NHM16` by entering explicit pieces of its first
visible exceptional divisor.  It proves nine exact generic or algebraic
slice exclusions.  It does not classify the complete divisor.

Replay all nine exact calculations with

```bash
.venv/bin/python scripts/verify_hc4_smooth_quartic_squarefree_line_exceptional_slices.py
```

or replay one slice with `--group` and one of the names printed by
`--help`.  SymPy constructs the same 81 reciprocal/Hessian coefficient
equations as `HC4NHM16`; Singular 4.4.1 computes the staged bases over
rational-function fields and, in three slices, exact quotient fields.  The
checker also recombines each linear and reduced basis and verifies zero
remainders for all 18 active variables (using $b_{12}^2$ on the one
set-theoretic slice).

## 1. Starting equations

Retain the notation of `HC4NHM16`:

\[
s_3=\frac{x^3+y^3}{3}+z^2(ux+vy)+wz^3,
\qquad
\ell=y+\tau x+\sigma z,
\tag{1.1}
\]

\[
A=A_0+zB,qquad
A_0=\begin{pmatrix}
0&0&-y^2\\
0&0&x^2\\
-y^2&x^2&p x^2+qxy+r y^2
\end{pmatrix}.
\tag{1.2}
\]

The 81 equations are the Hessian curls of

\[
C=\frac{\operatorname{adj}(A)+dd^{\mathsf T}}z,
\tag{1.3}
\]

the divisibility of \(R=\det(A)/z\) by \(\ell\), and

\[
\ell(R+e^{\mathsf T}d)=zR,
\qquad e=Ad/z.
\tag{1.4}
\]

At the generic point, the first linear reduction in `HC4NHM16` displays the
nonzero pivot

\[
\begin{aligned}
\Delta={}&(3p^2-qr)\tau^5+(9pr-q^2)\tau^4
 +(18r^2-6pq)\tau^3\\
&+(18p^2-6qr)\tau^2+(9pr-q^2)\tau+(3r^2-pq).
\end{aligned}
\tag{1.5}
\]

This note enters selected pieces of \(\Delta=0\).

## 2. The central quadratic vanishes

First set

\[
p=q=r=0.
\tag{2.1}
\]

Over
\(\mathbb Q(\tau,\sigma,b_{15},b_{16},b_{17})\), the pure linear basis has
eight pivots.  Reduction of the other equations gives

\[
(w,v,u,b_{14},b_{11},b_{10},b_9,b_8,b_5,b_2).
\tag{2.2}
\]

Together the two stages kill every active coefficient.  The first linear
basis visibly uses \(6\tau^3+1\).  On the algebraic subcase

\[
6\tau^3+1=0,
\tag{2.3}
\]

an exact quotient-ring computation gives instead

\[
(w,v,u,b_{14},b_{11},b_{10},b_9,b_8,b_5,b_4,b_2),
\tag{2.4}
\]

and the remaining linear pivots again kill every coefficient.  Thus this
first visible subexception of (2.1) is also empty.

## 3. The fiber \(\tau=0\)

Here

\[
\Delta=3r^2-pq.
\tag{3.1}
\]

On the chart \(p\ne0\), write

\[
p=3c,\qquad q=cm^2,\qquad r=cm,
\qquad c\ne0.
\tag{3.2}
\]

Over the total function field in
\(c,m,\sigma,b_{15},b_{16},b_{17}\), the reduced basis is

\[
(w,v,u,b_{14},b_{13},b_{11},b_8,b_5,b_2).
\tag{3.3}
\]

The nine linear pivots kill the complementary coefficients.  Its visible
secondary pivot is \(m^3-48\).  On the exact algebraic slice

\[
m^3=48,
\tag{3.4}
\]

the quotient-field linear basis first proves the nine substitutions used by
the checker.  Seven equations in the remaining core variables then have
reduced basis

\[
(v,u,b_{10}),
\tag{3.5}
\]

and, after this core reduction, the other equations have basis

\[
(w,b_{14},b_{11},b_8,b_5,b_2).
\tag{3.6}
\]

The checker reduces the nine substitution identities by the original linear
ideal before using them, then recombines all three stages and reduces the 18
original active coordinates to zero.  Thus the first visible secondary
slice is closed without the standard-basis timeout of the unstaged system.

On the chart \(p=0\), equation (3.1) gives \(r=0\).  With \(q\) generic,
the same reduced basis (3.3) occurs and the linear stage again kills the
complement.  The specialization \(q=0\) lies in Section 2.  Thus the generic
points of both standard charts of (3.1), as well as the first visible closed
subslice (3.4), are excluded.  Further factors not yet exposed by a complete
membership certificate may still define lower closed subsets.

## 4. The fiber \(\tau=-1\)

The first divisor splits:

\[
\Delta=5(p-r)(3p+q+3r).
\tag{4.1}
\]

### 4.1 The component \(p=r\)

At its generic point the reduced basis is

\[
(w,v,b_{14},b_{13},b_{11},b_8,b_5,b_2).
\tag{4.2}
\]

The linear stage contains \(u-v\) and pivots the other ten coefficients.
Its first visible secondary coefficient is

\[
E_1=q^2+3pq+8p^2.
\tag{4.3}
\]

On \(p\ne0\), put \(q=kp\) and impose

\[
k^2+3k+8=0.
\tag{4.4}
\]

The exact quotient-field reduced basis is

\[
(w,v,b_{14},b_{11},b_{10},b_9,b_8,b_5,b_2,b_{12}^2).
\tag{4.5}
\]

Hence $b_{12}=0$ set-theoretically in characteristic zero, and the linear
relations kill all other active coefficients.  The $p=0$ endpoint of
(4.3) lies in Section 2.

### 4.2 The component \(3p+q+3r=0\)

At its generic point the reduced basis is again (4.2).  The first visible
secondary pivot is

\[
E_2=7p-33r.
\tag{4.6}
\]

Substituting

\[
r=\frac7{33}p,qquad q=-\frac{40}{11}p
\tag{4.7}
\]

gives the exact reduced basis

\[
(w,v,b_{14},b_{12},b_{11},b_8,b_5,b_2).
\tag{4.8}
\]

The ten linear relations then kill the complement.  Thus both rational
components in (4.1), and their first displayed secondary-pivot strata, have
only determinant-zero support.

## 5. Result and remaining divisor

> **Theorem `HC4NHM17` -- Squarefree-line exceptional slices.**  In the
> squarefree-line row of `HC4NHM16`, the necessary reciprocal/Hessian system
> has only \(\det A=0\) support on each of the following nine strata:
> the generic locus $p=q=r=0$; its algebraic slice
> $6\tau^3+1=0$; the generic $p\ne0$ and $p=0$ charts of
> $\tau=0$, $3r^2-pq=0$; the algebraic slice $m^3=48$ on the
> $p\ne0$ chart; the two generic components of
> $\tau=-1$, $\Delta=0$; and the first visible secondary-pivot locus on each
> of those two components.  On the algebraic $E_1=0$ slice the certificate
> is set-theoretic because it contains $b_{12}^2$.

After every listed basis is combined with its linear stage, all
\(b_0,\ldots,b_{14},u,v,w\) vanish set-theoretically.  Therefore only the
bottom-right normal coefficients $b_{15},b_{16},b_{17}$ remain, and the
matrix has the determinant-zero form (3.7) of `HC4NHM16`.  No listed stratum
can carry a nonzero smooth quartic quotient.

What remains in this first approach is precise:

1. the generic divisor \(\Delta=0\) for arbitrary \(\tau\), not merely its
   fibers at \(0\) and \(-1\);
2. any additional factors hidden in nonlinear membership certificates,
   including further lower strata inside the displayed fibers;
3. the complementary residual-line chart and the other three
   basepoint-free boundary types.

The basepointed, doubled-line, and \(\mu=0\) rows remain separate.
