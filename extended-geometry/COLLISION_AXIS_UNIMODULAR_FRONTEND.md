# Collision-axis unimodular frontend

## Status

This note isolates the exact one-variable information carried by a normalized
three-dimensional Keller collision.  It proves a degree-independent pure-axis
support lower bound, gives the sharp example, and audits the proposed
elementary-length invariant.

The result gives the degree-independent global necessity of at least three
pure-axis nonlinear occurrences.  Its sharpness is only for admissible axis
data: no three-occurrence Keller collision is constructed.  It therefore does
not improve the stronger lower bound seven inside the existing degree-seven
census.  Completion over \(k[t]\) and first-jet integrability are automatic;
the first unresolved condition is extension to a three-variable Jacobian with
determinant one everywhere.

Throughout, \(k\) is a characteristic-zero field and integration from zero to
one means the coefficientwise polynomial antiderivative evaluated at those two
points.

## 1. Exact axis reduction

Let

\[
 F(0)=0,\qquad F(e_1)=0,\qquad JF(0)=I,\qquad \det JF=1,
\]

and put

\[
 h(t)=F(te_1),\qquad A(t)=JF(te_1),\qquad v(t)=h'(t)=A(t)e_1.
\]

Then

\[
 A(t)\in SL_3(k[t]),\qquad A(0)=I,
\]

and hence

\[
 v(0)=e_1,\qquad \int_0^1v(t)\,dt=h(1)-h(0)=0,
 \qquad (v_1,v_2,v_3)=k[t].
\tag{1.1}
\]

The last equality is the unimodularity condition.  It is stronger than the
collision equations alone and should be imposed before any off-axis Newton
support is introduced.

There is also an exact converse at the level of axis matrices.

> **Axis-completion lemma.**  Suppose \(v\in k[t]^3\) satisfies
> \(v(0)=e_1\).  There is an \(A\in SL_3(k[t])\) with \(A(0)=I\) and
> \(Ae_1=v\) if and only if \(v\) is unimodular.

Necessity is immediate.  For sufficiency, \(k[t]\) is a Euclidean domain, so
a unimodular column extends to a matrix \(B\in SL_3(k[t])\).  Since
\(B(0)e_1=e_1\), also \(B(0)^{-1}e_1=e_1\).  Therefore

\[
 A(t)=B(t)B(0)^{-1}
\]

has the required first column and satisfies \(A(0)=I\).  The same Euclidean
reduction shows \(SL_3(k[t])=E_3(k[t])\); an algorithmic reference for this
standard fact is Park and Woodburn,
[An Algorithmic Proof of Suslin's Stability Theorem over Polynomial Rings](https://arxiv.org/abs/alg-geom/9405003).

Thus no information beyond unimodularity is hidden in the existence of an
axis completion.

## 2. Why unrestricted elementary equivalence is too coarse

The action of \(E_3(k[t])\) is transitive on unimodular columns.  Consequently
there is only one unrestricted elementary orbit, not a useful classification
of collision columns.

Moreover, the constraints in (1.1) do not define an invariant subset for this
action.  Even the relative condition \(P(0)=I\) is insufficient: multiplication
by \(E_{21}(t)\) preserves the value at zero but generally changes the moment
\(\int_0^1v\).  Hence “classify (1.1) up to \(E_3(k[t])\)” is not a well-defined
quotient problem unless a substantially smaller moment-preserving groupoid is
specified.

The proposed word-length target also collapses.

> **Elementary-length proposition.**  With elementary length measured using
> factors \(E_{ij}(p(t))=I+p(t)e_{ij}\), the minimum length of an
> \(A\in SL_3(k[t])\) satisfying
> \[
> A(0)=I,\qquad \int_0^1A(t)e_1\,dt=0
> \]
> is exactly two.

A word of length one always has first entry equal to one in its first column,
so its first-column integral cannot vanish.  For the upper bound take

\[
 f(t)=24t,\qquad g(t)=t-\frac32t^2
\]

and

\[
 A(t)=E_{12}(f)E_{21}(g)
 =
 \begin{pmatrix}
 1+fg&f&0\\
 g&1&0\\
 0&0&1
 \end{pmatrix}.
\tag{2.1}
\]

Then \(A(0)=I\), \(\det A=1\), and

\[
 \int_0^1g(t)\,dt=0,
 \qquad
 \int_0^1(1+f(t)g(t))\,dt=1+8-9=0.
\]

Thus elementary length cannot by itself force a large mixed support in a
three-variable lift.

## 3. Sharp degree-independent axis-support bound

Count nonlinear pure-axis occurrences in \(h=F|_{ke_1}\), with coordinate
multiplicity.  Thus the fixed term \(t\) in \(h_1\) is not counted.

> **Pure-axis support theorem.**  If
> \[
> h(0)=h(1)=0,\qquad h'(0)=e_1,
> \qquad (h_1',h_2',h_3')=k[t],
> \]
> then \(h\) has at least three nonlinear monomial occurrences.  This bound is
> sharp in every characteristic-zero field.

Indeed, \(h_1(1)=0\) forces at least one nonlinear term in \(h_1\).  Neither
\(h_2\) nor \(h_3\) can have exactly one nonzero term, because its value at one
would then be nonzero.  With at most two nonlinear occurrences one must
therefore have \(h_2=h_3=0\).  Unimodularity would force \(h_1'\) to be a unit;
as \(h_1'(0)=1\), this gives \(h_1=t\), contradicting \(h_1(1)=0\).

Sharpness is witnessed by

\[
 h(t)=(t-t^2,\ t^2-t^3,\ 0),
 \qquad
 v(t)=(1-2t,\ 2t-3t^2,\ 0).
\tag{3.1}
\]

The two nonzero entries of \(v\) are coprime.  One normalized completion is

\[
 \begin{pmatrix}
 1-2t&-8t&0\\
 2t-3t^2&1+2t-12t^2&0\\
 0&0&1
 \end{pmatrix},
\tag{3.2}
\]

whose determinant is one and whose value at zero is the identity.
This witnesses sharpness of the one-variable frontend, not existence of a
three-variable Keller collision with three nonlinear terms.

Equality in the support bound has a simple classification, up to exchanging
the last two coordinates:

\[
 h(t)=\bigl(t-t^p,\ c(t^q-t^r),\ 0\bigr),
 \quad c\ne0,\quad p,q,r\ge2,\quad q<r.
\tag{3.3}
\]

The column in (3.3) is unimodular precisely when

\[
 \boxed{r^{p-1}\ne p^{r-q}q^{p-1}.}
\tag{3.4}
\]

To see this, a common nonzero root would satisfy

\[
 t^{p-1}=\frac1p,\qquad t^{r-q}=\frac qr.
\]

The standard compatibility condition for these two binomial equations is
exactly the equality excluded in (3.4).  The condition is substantive in
higher degree: \((p,q,r)=(9,2,6)\) fails because
\(6^8=9^4 2^8\).  Every triple with \(2\le p,q,r\le7\) and \(q<r\) passes.

The theorem gives a genuine degree-independent lower bound, but the sharp
value is three.  Any stronger total-support bound must use transverse/global
Keller compatibility, not the collision column alone.

## 4. First-jet integrability is automatic

Write a completed matrix as

\[
 A(t)=[v(t)\mid w(t)\mid u(t)]\in SL_3(k[t])
\]

and let \(H'(t)=v(t)\), \(H(0)=0\).  If \(v\) has zero integral, then
\(H(1)=0\).  The polynomial map

\[
 \widetilde F(t,y,z)=H(t)+y w(t)+z u(t)
\tag{4.1}
\]

satisfies

\[
 J\widetilde F(t,0,0)=A(t),\qquad
 \widetilde F(0)=\widetilde F(e_1)=0,
 \qquad J\widetilde F(0)=I.
\]

Thus there is no separate first-jet integrability obstruction.  However,

\[
 \det J\widetilde F
 =1+y\det[w',w,u]+z\det[u',w,u],
\tag{4.2}
\]

which need not equal one away from the axis.  Higher transverse terms must
cancel these defects and all later transverse coefficients while preserving
Hessian symmetry.  That global extension problem is the correct second
stage.

The resulting rigorous search order is therefore

\[
 \boxed{
 \text{axis restrictions }h
 \longrightarrow
 \gcd(h_1',h_2',h_3')=1
 \longrightarrow
 \text{normalized }SL_3(k[t])\text{ completion}
 \longrightarrow
 \text{global transverse Keller lift}
 \longrightarrow
 \text{remaining 3D support}.}
\tag{4.3}
\]

The completion arrow is constructive and carries no further obstruction.

## 5. Effect on the pinned degree-seven census

The support ledger in
[the global low-degree census](GLOBAL_LOW_DEGREE_SUPPORT_CENSUS.md) already
enforces the collision rule

\[
 a_1\ge1,\qquad a_2\ne1,\qquad a_3\ne1,
\]

where \(a_i\) counts nonlinear pure-axis terms in \(F_i\).  The unimodular
frontend adds

\[
 a_2+a_3>0.
\tag{5.1}
\]

Indeed, if \(a_2=a_3=0\), then
\(v=(h_1',0,0)\); zero integral prevents \(h_1'\) from being a unit.  Combined
with the collision rule, (5.1) says that at least one of \(a_2,a_3\) is at
least two.

Auditing the pinned determinant-balanced ledger gives:

| nonlinear support | balanced labelled supports | after axis gate | balanced orbits | after axis gate |
|---:|---:|---:|---:|---:|
| 4 | 30 | 0 | 15 | 0 |
| 5 | 85 | 0 | 47 | 0 |
| 6 | 1,694 | 900 | 851 | 450 |

Thus a reordered implementation can discard all size-four and size-five
balanced supports and nearly half of size six before coefficient algebra.
For the surviving minimum-axis patterns through degree seven, (3.4) is
automatic; in higher-degree searches it becomes an additional exact
resultant gate.

This pruning does not alter the proved lower bound seven: the original exact
coefficient computations already eliminated all 913 balanced orbit
representatives.  It improves the universal frontend and the routing of
future cardinality shards.

## 6. Reproduction

Run the exact symbolic audit with

```bash
.venv/bin/python scripts/verify_collision_axis_unimodular_frontend.py
```

The checker verifies (2.1), (3.1)--(3.4), the normalized completion (3.2),
the first-jet calculation (4.1)--(4.2), and the pinned census pruning counts.
