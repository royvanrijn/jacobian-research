# The squarefree quartic-denominator synchronization frontend

## Status

This note proves `HC4NHM9`, an exact linear frontend for the last clean
generic-corank-one quartic-denominator partition

\[
P=L_1L_2L_3L_4,
\qquad
\det\operatorname{Hess}(h_5)=P^2\ell,
\tag{0.1}
\]

where the \(L_i\) are distinct lines. It does **not** prove that the
partition is empty. It reduces that assertion to synchronization of four
constant polar directions, and separates the three line-arrangement types.

Replay the exact local identities with

~~~bash
.venv/bin/python scripts/verify_hc4_squarefree_quartic_denominator_frontend.py
~~~

The bounded finite-field sweep described in Section 4 is an experiment, not
a proof or a registered exclusion. The subsequent
[concurrence closure](HC4_SQUAREFREE_QUARTIC_CONCURRENCE_CLOSURE.md) proves
`HC4NHM10`: it closes the all-four-concurrent arrangement and the eight
exactly-three-concurrent patterns whose fourth flag is transverse. The
subsequent [general-position closure](HC4_SQUAREFREE_QUARTIC_GENERAL_POSITION_CLOSURE.md)
proves `HC4NHM11` and closes all sixteen no-three-concurrent patterns.

## 1. Pole order one forces a constant kernel

Let \(C=\operatorname{Hess}(h_5)\), let \(e=P C^{-1}d\) be the primitive
cleared vector, and fix an essential component \(L_i\). Its exponent in the
minimal denominator is

\[
b_i=1.
\]

The residue bound of `HC4NHM1` gives

\[
0\leq\kappa_i\leq b_i-1=0.
\tag{1.1}
\]

Thus the saturated kernel of \(C|_{L_i}\) is represented by a constant
nonzero vector \(v_i\). This conclusion uses neither `HC4NHM2` nor
`HC4NHM3`; it is special to the squarefree denominator.

## 2. Constant kernels are multiple-line polars

Write \(D_v\) for constant directional differentiation. Since

\[
C v_i=\nabla(D_{v_i}h_5),
\]

the constant-kernel condition \(C v_i\equiv0\pmod{L_i}\) is equivalent to

\[
\boxed{D_{v_i}h_5\in(L_i^2).}
\tag{2.1}
\]

Indeed, the restriction of the homogeneous quartic \(D_{v_i}h_5\) to the
line has zero gradient, hence is zero; its first normal derivative is zero
as well. Conversely, (2.1) makes every entry of \(Cv_i\) divisible by
\(L_i\).

There are two local rows.

### 2.1 Tangent kernel

If \(L_i(v_i)=0\), normalize \(L_i=x\) and \(v_i=\partial_z\). Equation
(2.1) gives

\[
h_5=F_5(x,y)+x^2H_3(x,y,z).
\tag{2.2}
\]

The last Hessian column has orders \((x,x^2,x^2)\). Symmetry then gives

\[
x^2\mid\det C
\tag{2.3}
\]

automatically. Thus a tangent constant kernel needs no further local
determinant equation.

### 2.2 Transverse kernel

If \(L_i(v_i)\ne0\), normalize \(L_i=x\) and
\(v_i=\partial_x\). Equation (2.1) gives

\[
h_5=F_5(y,z)+x^3J_2(y,z)+x^4K_1(y,z)+c x^5.
\tag{2.4}
\]

Its first determinant coefficient is

\[
[x]\det C
=6J_2\det\operatorname{Hess}_{y,z}(F_5).
\tag{2.5}
\]

On the generic-corank-one locus the binary Hessian in (2.5) is nonzero.
Consequently

\[
x^2\mid\det C
\quad\Longleftrightarrow\quad
J_2=0
\quad\Longleftrightarrow\quad
D_{v_i}h_5\in(x^3).
\tag{2.6}

The binary-Hessian-zero alternative is a lower-Smith boundary and is not
absorbed into this theorem.

Combining the two rows, every packet in (0.1) satisfies the four **linear**
polar conditions

\[
D_{v_i}h_5\in
\begin{cases}
(L_i^2),&L_i(v_i)=0,\\
(L_i^3),&L_i(v_i)\ne0.
\end{cases}
\tag{2.7}

Once the flags \((L_i,v_i)\) are fixed, no determinant expansion is needed
to enforce \(P^2\mid\det C\). Because both sides have degrees eight and
nine, respectively, every nonzero survivor automatically has
\(\det C=P^2\ell\) for a line \(\ell\).

## 3. Finite geometric partition

Four distinct projective lines have exactly three incidence types:

1. no three concurrent, normalized to
   \((x,y,z,x+y+z)\);
2. exactly three concurrent, normalized to \((x,y,x+y,z)\);
3. all four concurrent, represented by a binary squarefree quartic line
   pencil \((x,y,x+y,x+\lambda y)\), with the usual cross-ratio parameter;
   \(\lambda=-1\) is one convenient chart.

For each arrangement there are sixteen tangent/transverse flag patterns.
Within a pattern, (2.7) is a matrix of linear equations in the twenty-one
coefficients of \(h_5\), with at most eight projective kernel parameters.
This is the correct next elimination: classify the parameter values at which
that matrix has a non-cone nullspace, then impose the cleared module and
gradient equations. It is much smaller than eliminating the coefficients of
an arbitrary ternary quintic Hessian determinant.

There is also an immediate collision rule for repeated kernel directions.
If \(v_i=v_j=v\), then \(D_vh_5\) is divisible by both prescribed powers of
the distinct lines \(L_i,L_j\). A tangent--transverse or
transverse--transverse pair has total required degree greater than four, so
\(D_vh_5=0\) and the Hessian determinant vanishes. A nonzero determinant can
therefore have a repeated direction only on a tangent--tangent pair, where

\[
D_vh_5=cL_i^2L_j^2.
\tag{3.1}

No direction can occur on three denominator lines.

## 4. Bounded synchronization experiment

The exploratory script

~~~bash
.venv/bin/python scripts/research_hc4_squarefree_quartic_kernel_sync.py \
  --arrangement general --finite-field-prime 3
.venv/bin/python scripts/research_hc4_squarefree_quartic_kernel_sync.py \
  --arrangement triple --finite-field-prime 3
.venv/bin/python scripts/research_hc4_squarefree_quartic_kernel_sync.py \
  --arrangement pencil --finite-field-prime 3
~~~

exhausts all \(13^4=28,561\) quadruples of projective kernel directions over
\(\mathbf F_3\) for each displayed arrangement. Many synchronized linear
spaces survive. Testing three deterministic nullspace assignments at three
source points in every such space finds no nonzero Hessian determinant
witness. The pencil command covers the displayed \(\lambda=-1\) chart.

This is evidence for the prospective synchronization theorem

\[
\text{four conditions (2.7)}\quad\Longrightarrow\quad\det C=0,
\tag{4.1}
\]

not a proof of it in characteristic zero or even an exhaustive coefficient
census over \(\mathbf F_3\). In particular, an algebraic kernel direction may
be defined only over an extension field, the pencil cross ratio varies, and
deterministic witness testing does not establish ideal containment over
\(\mathbf Q\).

> **Theorem `HC4NHM9` -- Squarefree synchronization frontend.** On the clean
> generic-corank-one packet with squarefree quartic minimal denominator, all
> four saturated kernel directions are constant and satisfy the linear polar
> divisibilities (2.7). The packet splits into the three line arrangements
> and sixteen tangent/transverse patterns above. The squarefree partition is
> reduced to this flag synchronization problem.

The subsequent theorems `HC4NHM10--12` prove (4.1) on all forty-eight rows.
The last tangent-fourth calculation is
[`HC4_SQUAREFREE_QUARTIC_TANGENT_FOURTH_CLOSURE.md`](HC4_SQUAREFREE_QUARTIC_TANGENT_FOURTH_CLOSURE.md).
