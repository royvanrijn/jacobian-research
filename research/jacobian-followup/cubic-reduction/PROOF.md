# A cubic Keller collision in seventeen variables

**Historical dominated example.** Twelve-variable cubic maps had already
been verified here before this experiment. The [baseline correction](BASELINE_CORRECTION.md)
records the missed comparison. The current support follow-up is the
[36-term twelve-variable conjugate](PROOF12.md). This existence theorem
and its exact evidence are retained; it is not a dimension improvement.

There is a rational polynomial map \(H:\mathbb A^{17}\to\mathbb A^{17}\)
of ordinary total degree three with

\[
H(0)=0,\qquad JH(0)=I,\qquad \det JH=1,
\]

and three distinct rational points in one fibre. It has **58 nonlinear
monomial occurrences** in the displayed coordinates.

This removes two redundant carrier directions from the
[19-variable construction posted by rhicksrad, 20 July 2026](https://rhicksrad.github.io/jacobian-degree3/).
The foundational counterexample is Alpöge's; the stable cubic reduction
is the classical Bass–Connell–Wright method. No new counterexample
mechanism, priority claim, or optimal dimension is asserted.

Together with the quadratic midpoint identity below, this gives exact
minimum ordinary degree **three for characteristic-zero Keller
collisions when dimension may vary**. The minimum dimension for a cubic
example remains open. This does not settle degree minimization in three
variables or the three-variable normalized support minimum.

## Explicit map

Use the variables

\[
(x,y,z,w_1,\ldots,w_{10},u,w_{13},v,w_{16}),
\quad L=xy+w_9,\quad M=x^2+w_{13}.
\]

The first three outputs are

\[
\begin{aligned}
H_1={}&x-\tfrac32x^2y+\tfrac12Mw_{16}
                  +\tfrac12xz\,w_{13},\\
H_2={}&y+3xz+12xy^2+3Mv-6Mw_4-3w_2(w_4+yz)\\
 &+3w_5(xw_4-3w_7+w_8-3y^2)
   +3yw_{13}(3w_7-w_8-2z),\\
H_3={}&z+4y^2+3xyz+L(u-7w_7+3w_8)\\
 &+w_9(3xw_4-xw_6+6yw_3-7y^2)
   +(w_7+y^2)w_{10}+w_2(zw_7-w_1)\\
 &+w_5(xw_1-3yz-3w_4+w_6)-3w_3^2.
\end{aligned}
\]

The fourteen carrier outputs, in the same order as their source variables,
are

\[
\begin{aligned}
C_1&=w_1+y^2z,&
C_2&=w_2-xyw_{13}-x^2w_9-w_9w_{13},\\
C_3&=w_3+xy^2,&C_4&=w_4+yz,\\
C_5&=w_5+x^2y,&C_6&=w_6+xw_1,\\
C_7&=w_7+y^2,&C_8&=w_8+xw_4,\\
C_9&=w_9+xy,&C_{10}&=w_{10}+zw_2,\\
C_u&=u+6yw_3-xw_6,&C_{13}&=w_{13}+x^2,\\
C_v&=v+3yw_7-yw_8,&C_{16}&=w_{16}+xz.
\end{aligned}
\]

The nonlinear row counts are
\((4,15,21,1,3,1,1,1,1,1,1,1,1,2,1,2,1)\).
Their sum is 58. Degree three is visible in the formulas.

## Why two carriers can be removed

In the external 19-variable map \(G\), its first three outputs depend on
\(w_{11},w_{12}\) only through \(6w_{11}-w_{12}\), and on
\(w_{14},w_{15}\) only through \(3w_{14}-w_{15}\). No other retained
output uses these four variables except their own carrier rows.

Set

\[
u=6w_{11}-w_{12},\quad r=w_{11},\qquad
v=3w_{14}-w_{15},\quad s=w_{14}.
\]

Apply the same linear coordinate changes in the carrier targets:
retain \(6C_{11}-C_{12}\) and \(3C_{14}-C_{15}\), and retain
\(C_{11}\) and \(C_{14}\) as complementary outputs. The transformed
map has the form

\[
\bigl(\widetilde H(X),\ r+yw_3,\ s+yw_7\bigr),
\]

where all seventeen coordinates of \(\widetilde H\) are independent of
\(r,s\). Source shears turn the last two outputs into identity
coordinates without changing \(\widetilde H\). Thus \(G\) is stably
equivalent to \((\widetilde H,\mathrm{id}_{\mathbb A^2})\).
Finally replace its first three outputs by
\((\widetilde H_3/2,\widetilde H_2,\widetilde H_1)\) to get the displayed
map \(H\). This gives the identity derivative and determinant one.

The independent checker verifies all nineteen polynomial substitution
identities, the retained rows' independence of the complementary
directions, and both removed output formulas. It does not assume that
the external construction has been reviewed or correctly verified.

## Full determinant and collision checks

Write the seventeen-by-seventeen Jacobian in blocks with the original
three coordinates first:

\[
JH=\begin{pmatrix}A&B\\C&D\end{pmatrix},\qquad D=I+N.
\]

The complete polynomial calculation gives \(N^3=0\). Therefore
\(D^{-1}=I-N+N^2\), whose two inverse identities are checked exactly.
The checker forms the complete polynomial Schur complement
\(A-BD^{-1}C\) and its three-by-three determinant. It is exactly one.
All variables, including every carrier, remain symbolic in this check.

[The coefficient witness](cubic17-witness.json) contains three rational
points, beginning with source coordinates
\((0,0,-1/4)\), \((1,-3/2,13/2)\), and \((-1,3/2,13/2)\).
The full points are distinct and all map to
\((0,0,-1/4,0,\ldots,0)\). The checker evaluates every coordinate exactly.
The common fibre is not the origin fibre; the 58-term count refers to the
displayed tangent-normalized coordinates. An affine collision normalization
preserves degree three but need not preserve that term count.

## The degree lower bound in every dimension

For any polynomial map \(Q\) of degree at most two over a field of
characteristic zero, polarization gives

\[
Q(a)-Q(b)=JQ((a+b)/2)(a-b).
\]

A nonzero constant Jacobian determinant makes the matrix on the right
invertible. Thus \(Q(a)=Q(b)\) forces \(a=b\). The identity applies to
each constant, linear, square and mixed quadratic monomial, in any
dimension. The independent checker includes the square and mixed
quadratic identities in rational sparse arithmetic.

No degree-at-most-two Keller collision exists, whereas the displayed
degree-three map has one. This establishes the stated exact degree
threshold, independently of any dimension minimum.

## Reproduction and search boundaries

~~~sh
.venv/bin/python research/jacobian-followup/cubic-reduction/construct17.py
.venv/bin/python research/jacobian-followup/cubic-reduction/verify.py
.venv/bin/python research/jacobian-followup/cubic-reduction/check_integration.py
~~~

The [primary construction](construct17.py) uses SymPy 1.14.0.
The [independent checker](verify.py) uses only Python's standard library,
with rational sparse polynomial arithmetic. The
[validation receipt](validation.json) pins all inputs.
[Provenance](sources/provenance.json) pins the external machine-readable
map by URL and SHA256.

Earlier bounded discovery attempts are retained in
[the exploration record](EXPLORATION.md). None supplies a lower bound
on dimension or term count. The normalized three-variable degree-seven
support interval remains **[10,48]**. The global minimum dimension,
the cubic minimum dimension, and support optimality remain open.

