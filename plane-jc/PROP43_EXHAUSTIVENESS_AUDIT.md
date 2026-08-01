# Finite exhaustiveness audit for GGHV Proposition 4.3

## Result and claim boundary

This note fills the two finite arithmetic steps abbreviated in Proposition 4.3
of Guccione--Guccione--Horruitiner--Valqui (GGHV),
[*Increasing the degree of a possible counterexample to the Jacobian Conjecture
from 100 to 108*](https://arxiv.org/abs/2204.14178):

1. the phrase “as in Proposition 4.1”, which passes from the reduced corner
   `(24,7)` to the four listed alternatives
   `(24,7),(17,5),(10,3),(3,1)` and successor direction `(-2,7)`; and
2. the sentence excluding `k=2` because the two terminal edges “would have no
   way of being parallel”.

The audit proves, by exact integer arithmetic, that the cited general
Newton-polygon results force precisely the two corrected Laurent polygons used
by the `(72,108)` certificate replay.  It also identifies the apparent
`(0,1)`/`(8,14)` typo in the displayed statement of Proposition 4.3.

The executable checker is
[`cas/verify_prop43_exhaustiveness.py`](cas/verify_prop43_exhaustiveness.py)
and uses only the Python standard library.

This is not a new proof of the general GGV/GGHV structural theorems or of the
earlier admissible-chain census.  Those remain cited theorem inputs.  The new
content is the complete case-specific bridge from their outputs to the two
Laurent systems.

## 1. Published theorem interface

The proof starts from the following published outputs.

- The last chain is
  `A0=(8,28), A1=(11/4,7), (m,n)=(3,2)`.
- After the coordinate swap and Corollary 7.4, the normalized upper corner is
  `(28,8)`, the root degree is `q=4`, and the auxiliary endpoint is `(7,2)`.
- Proposition 3.12 and Proposition 2.5 leave predecessor directions
  `(1,-2)` and `(1,-3)`.
- The displayed Laurent triangular transformations leave the three normalized
  intermediate polygons

  ```text
  a) {(-2,0),(0,0),(28,8),(0,1)}
  b) {(-3,0),(0,0),(28,8),(0,1)}
  c) {(-3,0),(0,0),(16,4),(28,8),(0,1)}.
  ```

- Proposition 8.2 of Guccione--Guccione--Valqui,
  [*On the shape of possible counterexamples to the Jacobian Conjecture*](https://arxiv.org/abs/1401.1784),
  supplies the dichotomy used below.  If the next endpoints remain
  proportional, the normalized endpoint `(a,b)` is integral, has smaller old
  weight, and has positive determinant with the current corner.  At the first
  nonproportional successor there is a positive integer `k` with

  ```text
  (k+1)b < a,
  {en(P),en(Q)} = {(-k,0),(k+1,1)}.
  ```

The small predecessor calculation can also be replayed directly.  Its integer
parameter satisfies

\[
1<\Delta<\frac72,
\qquad
7-2\Delta\mid\Delta-1.
\]

Only `Delta=3` survives, giving `(1,-3)` beside the separately isolated
exceptional direction `(1,-2)`.

## 2. Exact upper-edge reduction

The common upper edge is represented by

\[
y(x^4y-\alpha)^7.
\]

Under `y -> y + alpha*x^-4`,

\[
\begin{aligned}
(y+\alpha x^{-4})
\bigl(x^4(y+\alpha x^{-4})-\alpha\bigr)^7
&=(y+\alpha x^{-4})(x^4y)^7\\
&=x^{28}y^8+\alpha x^{24}y^7.
\end{aligned}
\]

Thus the edge `(28,8)--(0,1)` is reduced exactly to
`(28,8)--(24,7)`.  There is no genericity assumption or unexamined
cancellation in this step.

## 3. A uniform prime-defect lemma

The finite calculation in Proposition 4.1 has a useful uniform form.
Let the current normalized corner be `(A,B)`, and let a proper proportional
endpoint be `(a,b)`.  Put

\[
g=\gcd(A-a,B-b),
\qquad
(p,q)=\frac1g(A-a,B-b).
\]

The auxiliary homogeneous element has its second corner at

\[
(1,1)+c(p,q),
\qquad c\in\mathbb N,
\]

and this point lies on the ray through `(A,B)`.  Applying the normal
`(-B,A)` gives the necessary divisibility

\[
Ab-Ba\mid(A-B)g.
\]

Since

\[
Ab-Ba=g(Bp-Aq),
\]

we obtain the scale-free condition

\[
\boxed{Bp-Aq\mid A-B.}
\]

For the printed Proposition-4.1 corner `(A,B)=(21,8)`, the defect is the prime
`13`.  The admissible primitive step of defect one is `(8,3)`, giving the two
proper scales `(13,5)` and `(5,2)`; defect thirteen gives the forbidden
diagonal point `(1,1)`.  This reproduces the published table.

For the live corner `(A,B)=(24,7)`, the defect is the prime `17`.  The two
possible primitive equations are

\[
7p-24q=1
\quad\text{or}\quad
7p-24q=17.
\]

Inside the displayed polygon envelope they have respectively

\[
(p,q)=(7,2)
\quad\text{and}\quad
(p,q)=(23,6).
\]

The first direction permits scales `g=1,2,3`, giving

\[
(17,5),\qquad(10,3),\qquad(3,1),
\]

whereas the second gives only the forbidden diagonal point `(1,1)`.
The checker also performs the complete lattice census, described next, so the
Diophantine simplification is not being used as an unverified shortcut.

## 4. Complete opposite-corner census

Let the current corner be `A=(24,7)` and the old direction be `(-1,4)`.
Every proportional next endpoint in any of the three displayed intermediate
polygons lies in the safe common envelope

\[
-3\le a\le24,
\qquad
0\le b\le7,
\]

and Proposition 8.2(1) requires

\[
v_{-1,4}(a,b)<v_{-1,4}(24,7)=4,
\qquad
24b-7a>0.
\]

Exact enumeration gives twelve points:

\[
\begin{split}
&(-3,0),(-2,0),(-1,0),(1,1),(2,1),(3,1),\\
&(5,2),(6,2),(9,3),(10,3),(13,4),(17,5).
\end{split}
\]

The auxiliary-element condition is

\[
\boxed{24b-7a\mid17\gcd(24-a,7-b).}
\]

The complete table is:

| `(a,b)` | `24b-7a` | `gcd(24-a,7-b)` | `17g` | result |
|---:|---:|---:|---:|:---|
| `(-3,0)` | 21 | 1 | 17 | excluded |
| `(-2,0)` | 14 | 1 | 17 | excluded |
| `(-1,0)` | 7 | 1 | 17 | excluded |
| `(1,1)` | 17 | 1 | 17 | diagonal, forbidden |
| `(2,1)` | 10 | 2 | 34 | excluded |
| `(3,1)` | 3 | 3 | 51 | survives |
| `(5,2)` | 13 | 1 | 17 | excluded |
| `(6,2)` | 6 | 1 | 17 | excluded |
| `(9,3)` | 9 | 1 | 17 | excluded |
| `(10,3)` | 2 | 2 | 34 | survives |
| `(13,4)` | 5 | 1 | 17 | excluded |
| `(17,5)` | 1 | 1 | 17 | survives |

Thus the proper proportional alternatives are exactly

\[
(17,5),\qquad(10,3),\qquad(3,1).
\]

Together with the immediate nonproportional break at `(24,7)`, this reproduces
exactly the four alternatives printed in Proposition 4.3.

## 5. The common successor direction

All four alternatives satisfy

\[
\boxed{2a=7b-1},
\qquad\text{equivalently}\qquad
v_{-2,7}(a,b)=1.
\]

Moreover,

\[
(24,7)-(17,5)=(7,2),
\]

\[
(24,7)-(10,3)=2(7,2),
\]

\[
(24,7)-(3,1)=3(7,2).
\]

Hence every surviving proper continuation lies on the primitive edge direction
`(7,2)`, with primitive normal `(-2,7)`.  This proves the case-specific claim

\[
\operatorname{Succ}_P(-1,4)=
\operatorname{Succ}_Q(-1,4)=(-2,7)
\]

unless the nonproportional break occurs immediately at `(24,7)`; the latter
case has the same terminal direction by the calculation below.

## 6. Exact parallelism calculation

At a normalized corner `(a,b)`, the starts of the two terminal edges are
`2(a,b)` for `P` and `3(a,b)` for `Q`.  Proposition 8.2 gives the endpoint set

\[
\{(-k,0),(k+1,1)\}.
\]

For the assignment

\[
\operatorname{en}(P)=(-k,0),
\qquad
\operatorname{en}(Q)=(k+1,1),
\]

the parallelism determinant is

\[
D_-=-2a+(5k+2)b-k.
\]

For the swapped assignment it is

\[
D_+=3a-(5k+3)b+k.
\]

Using `2a=7b-1`, these simplify uniformly to

\[
\boxed{D_-=(k-1)(5b-1)},
\]

\[
\boxed{D_+=\frac{(3-2k)(5b-1)}2}.
\]

Since `b>=1`, `5b-1` is nonzero.  Therefore:

- for `k=1`, exactly the first assignment is parallel;
- for `k=2`, neither assignment is parallel.

Thus `k=2` is impossible, and the unique terminal assignment is

\[
\boxed{
\operatorname{en}(P)=(-1,0),
\qquad
\operatorname{en}(Q)=(2,1).}
\]

For `(3,1)`, `k=2` is already ruled out by `(k+1)b<a`; the determinant formula
agrees with that restriction.

The exact numerical table is:

| `(a,b)` | `k` | `D_-` | `D_+` |
|---:|---:|---:|---:|
| `(24,7)` | 1 | 0 | 17 |
| `(24,7)` | 2 | 34 | -17 |
| `(17,5)` | 1 | 0 | 12 |
| `(17,5)` | 2 | 24 | -12 |
| `(10,3)` | 1 | 0 | 7 |
| `(10,3)` | 2 | 14 | -7 |
| `(3,1)` | 1 | 0 | 2 |

## 7. Why the proper proportional alternatives are contradictions

For `k=1`, the `P` edge leaving `2(a,b)` is

\[
(-1,0)-2(a,b)=(-7b,-2b)=-b(7,2).
\]

The `Q` edge leaving `3(a,b)` is

\[
(2,1)-3(a,b)
=-\frac{3b-1}{2}(7,2).
\]

For each proper candidate `(17,5),(10,3),(3,1)`, the incoming edge from
`(24,7)` is already a negative multiple of `(7,2)`.  The alleged endpoint
`2(a,b)` or `3(a,b)` therefore lies strictly inside the same straight Newton
edge: the incoming and outgoing segments have the same ray and cannot define
two successive boundary directions.  This contradicts the assumption that
`(a,b)` was the endpoint of a distinct proportional successor edge.

Consequently the sole genuine nonproportional break occurs directly at
`(24,7)`, with successor normal `(-2,7)` and terminal endpoints `(-1,0)` for
`P` and `(2,1)` for `Q`.

## 8. Final monomial map and the two polygons

The final map is

\[
x\longmapsto x^{-1},
\qquad
y\longmapsto x^4y,
\]

and sends an exponent `(i,j)` to

\[
T(i,j)=(-i+4j,j).
\]

For cases `a)` and `b)`, the pre-map vertices are

\[
\begin{aligned}
N(P)&=\{(-1,0),(0,0),(56,16),(48,14)\},\\
N(Q)&=\{(2,1),(0,0),(84,24),(72,21)\}.
\end{aligned}
\]

Their images are

\[
\begin{aligned}
N(P)&=\{(0,0),(1,0),(8,14),(8,16)\},\\
N(Q)&=\{(0,0),(2,1),(12,21),(12,24)\}.
\end{aligned}
\]

Case `c)` additionally contains `2(16,4)=(32,8)` and
`3(16,4)=(48,12)`, which map to `(0,8)` and `(0,12)`.  It therefore gives

\[
\begin{aligned}
N(P)&=\{(0,0),(1,0),(8,14),(8,16),(0,8)\},\\
N(Q)&=\{(0,0),(2,1),(12,21),(12,24),(0,12)\}.
\end{aligned}
\]

The coordinate Jacobian of the monomial map is `-x^2`.  The chain rule changes
the preceding constant bracket into a nonzero scalar multiple of `x^2`, and a
scalar normalization gives `[P,Q]=x^2`.

## 9. Printed typo in Proposition 4.3

The displayed Case-1 statement in the arXiv PDF lists `(0,1)` where its own
final computation lists `(8,14)`.  The exponent map above independently forces
`(8,14)`.  In addition, `(0,1)` lies in the interior of the vertical segment
from `(0,0)` to `(0,8)`, so it cannot replace the missing lower-right vertex.
The corrected Case-1 polygon is therefore

\[
N(P)=\{(0,0),(1,0),(8,14),(8,16),(0,8)\}.
\]

The repository and exact-certificate archive already use this corrected form.

## 10. Consequence for the `(72,108)` programme

Assuming the cited general GGV/GGHV structural results and the previously
audited admissible chain, the exact finite calculations above force precisely
the two corrected Proposition-4.3 Laurent systems.  There is no omitted
opposite-corner branch, no `k=2` escape, and no third polygon.

The repository has separately replayed characteristic-zero unit certificates
excluding both Laurent systems.  Combining that replay with this finite
front-end audit closes the previously recorded interface gap in the
`(72,108)` exclusion.  The resulting degree-125 statement still inherits the
ordinary mathematical dependency on the cited theorem ladder and its earlier
bounded chain enumeration; it is not a from-scratch reproof of those papers or
of the full plane Jacobian conjecture.

## Reproduction

Run from the repository root:

```bash
python3 plane-jc/cas/verify_prop43_exhaustiveness.py
```

Expected final marker:

```text
PROP43_EXHAUSTIVENESS_AUDIT_PASS
```
