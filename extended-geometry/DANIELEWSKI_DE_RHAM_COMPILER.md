# The Danielewski de Rham compiler

The coupled volume equation has a two-dimensional cohomological shadow
that can be computed before any coefficient-ideal elimination.  This note
constructs that compiler explicitly for

\[
D=\{cw=P(x)\},\qquad P=x^3+x.
\]

It turns every polynomial two-form \(h\omega\) into two scalar coordinates
and reduces the maximal quadratic flux ledger to five equations.

## 1. Poisson-homology presentation

Let

\[
B=k[c,x,w]/(cw-P(x)),\qquad
\omega=\frac{dc\wedge dx}{c}.
\]

Since

\[
dF\wedge dG=\{F,G\}\omega,
\]

the top algebraic de Rham group is

\[
H^2_{\rm dR}(D)
\simeq
B/\operatorname{span}_k\{\{F,G\}:F,G\in B\}.        \tag{1}
\]

Every element of \(B\) has a normal form

\[
h_0(x)
+\sum_{i\ge1}c^i h_i(x)
+\sum_{i\ge1}w^i\widetilde h_i(x).                 \tag{2}
\]

The terms involving \(c\) or \(w\) vanish in (1), because

\[
\{x,c^i h(x)\}=-i c^i h(x),\qquad
\{x,w^i h(x)\}=i w^i h(x).                         \tag{3}
\]

For the remaining polynomial part,

\[
\{c,wf(x)\}
=P'(x)f(x)+P(x)f'(x)
=\bigl(P(x)f(x)\bigr)'.                            \tag{4}
\]

Conversely, reducing a general bracket to the normal form (2) shows that
its pure-\(x\) part lies in \(\partial_x(Pk[x])\).  Hence

\[
\boxed{
H^2_{\rm dR}(D)
\simeq
k[x]/\partial_x(Pk[x])
}                                                    \tag{5}
\]

as vector spaces.

For a squarefree polynomial of degree \(n\), the right side has basis

\[
1,x,\ldots,x^{n-2}.
\]

This recovers \(\dim H^2_{\rm dR}(D)=n-1\).

## 2. Cubic recurrence

For \(P=x^3+x\),

\[
\bigl(Px^m\bigr)'
=(m+3)x^{m+2}+(m+1)x^m.
\]

Therefore

\[
\boxed{
[x^{m+2}]
=-\frac{m+1}{m+3}[x^m].
}                                                    \tag{6}
\]

Every class has a unique form

\[
\rho(h)=\rho_0(h)[1]+\rho_1(h)[x].                 \tag{7}
\]

For example,

\[
\begin{array}{c|cccccc}
h&1&x&x^2&x^3&x^4&x^5\\ \hline
\rho(h)&1&x&-\frac13&-\frac12x&\frac15&\frac13x.
\end{array}
\]

The relation

\[
\rho(P')=\rho(3x^2+1)=0
\]

is the first consistency check.

## 3. Algorithm

To compile \(h\in B\):

1. replace every factor \(c^iw^j\) by
   \(c^{i-m}w^{j-m}P(x)^m\), where \(m=\min(i,j)\);
2. discard every remaining monomial containing \(c\) or \(w\), by (3);
3. reduce each pure power of \(x\) recursively with (6).

The result is the pair

\[
\boxed{\rho(h)=(\rho_0(h),\rho_1(h)).}              \tag{8}
\]

The procedure is linear, exact over \(\mathbb Q\), and avoids integration
or numerical periods.

## 4. Compiled coupled flux

For three target coordinates \(A,F,G\in B[a]\), put

\[
\eta=A\,d_DF\wedge d_DG
=A\{F,G\}\omega.
\]

The flux identity gives

\[
\frac d{da}\rho\!\left(A\{F,G\}\right)=(1,0)        \tag{9}
\]

whenever

\[
A_a\{F,G\}-F_a\{A,G\}+G_a\{A,F\}=1.
\]

Equivalently,

\[
\rho\!\left(A\{F,G\}\right)
=(a+\lambda,\mu)                                   \tag{10}
\]

for constants \(\lambda,\mu\).

Apply this to the normalized maximal ambient-quadratic ansatz

\[
(A,F,G)=(a,c,w)
+\text{linear combinations of }
(ac,ax,aw,a^2,c^2,x^2,w^2,cw,cx,xw).
\]

Although the pointwise determinant has 102 source-monomial coefficients,
the compiled flux has degrees

\[
\deg_a\rho_0\le2,\qquad
\deg_a\rho_1\le3.
\]

After ignoring the two arbitrary constant terms, (10) gives exactly

\[
\boxed{5\text{ nonconstant flux equations}.}       \tag{11}
\]

They consist of the coefficients

\[
[\rho_0]_{a},\ [\rho_0]_{a^2},\
[\rho_1]_{a},\ [\rho_1]_{a^2},\ [\rho_1]_{a^3},
\]

with the first equal to one and the other four equal to zero.

These five equations are exact linear combinations, in cohomology, of the
full determinant ledger.  Adding them as redundant generators did not make
the original row-ordered rational Groebner basis terminate within the
exploratory time window.  A separate column-major ordering subsequently
gave the exact rational unit basis and
[closed the full quadratic ansatz](DANIELEWSKI_FULL_QUADRATIC_OBSTRUCTION.md).
The compiler's continuing value is structural: every higher-degree ansatz
receives sparse flux gates before pointwise elimination.

## 5. Next use

The compiler enables two concrete continuations.

First, apply (8) to each coefficient in the leading-\(a\)-face descent.
The highest Poisson-rank-one face must have zero flux, while the next face
must contribute exactly \((1,0)\).  This is the cleanest route to an
all-degree contradiction.

Second, use the same compiled equations as prefilters for ambient degree
three and above.  The quadratic ideal no longer needs branching: its exact
rational unit basis is known.

## 6. Reproduction

Run

```bash
.venv/bin/python scripts/verify_danielewski_de_rham_compiler.py
```

The checker verifies the normal-form reduction, the recurrence, vanishing
on exact brackets, independence of the two period coordinates, the flux
identity, and the five-equation maximal quadratic compilation.
