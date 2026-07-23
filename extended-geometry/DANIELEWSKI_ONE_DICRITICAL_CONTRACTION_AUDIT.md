# One-dicritical Danielewski contraction audit

The first smooth Danielewski target with a nonempty missing boundary is

\[
X=\mathbb A^1_a\times D,\qquad
D=\{cw=x(x+1)\}.
\]

It contains affine three-space through

\[
\Phi(a,b,c)=(a,c,bc,b(bc+1)).
\]

The complement of the image is the single plane

\[
E=\{c=0,\ x=-1\}\simeq\mathbb A^2_{a,w}.
\]

The evident contraction

\[
(a,c,x)\colon X\longrightarrow\mathbb A^3
\]

contracts \(E\), but its pulled-back ordinary volume form is

\[
da\wedge dc\wedge dx=c\,\Omega_X,
\qquad
\Omega_X=\frac{da\wedge dc\wedge dx}{c}.
\]

Thus it has a Jacobian divisor.  The Keller-contraction problem asks instead
for \(A,F,G\in k[X]\) satisfying

\[
dA\wedge dF\wedge dG=\lambda\Omega_X,
\qquad \lambda\in k^\times.                       \tag{1}
\]

This note closes the complete ambient-quadratic layer of (1) and checks that
the announced foundational Keller map does not already solve (1) through a
mere permutation of its source coordinates.

## 1. Complete ambient-quadratic obstruction

Give \(a,c,x,w\) ambient degree one.  Modulo

\[
cw=x^2+x,
\]

the degree-at-most-two coordinate space has basis

\[
\begin{aligned}
1,\ a,\ c,\ x,\ w,\quad
&a^2,c^2,x^2,w^2,\\
&ac,ax,aw,cx,xw.                                  \tag{2}
\end{aligned}
\]

### Theorem 1.1

Over a field of characteristic zero, there are no three elements of the
space (2) satisfying (1).

### Proof

Pull (2) back by \(\Phi\), so

\[
x=bc,\qquad w=b(bc+1).
\]

If (1) holds, the ordinary source Jacobian is a nonzero constant.  Its
linear jet at the origin is therefore invertible.  The only basis elements
in (2) with nonzero source-linear terms are \(a,c,w\).  After a target
translation, a constant target linear change, and rescaling the desired
determinant, every candidate has the unique normalized shape

\[
(A,F,G)=(a,c,w)+
v_xx+\sum_N v_NN,                                  \tag{3}
\]

where \(N\) runs over the ten quadratic monomials in (2) and each coefficient
vector lies in \(k^3\).  Thus (3) has 30 scalar coefficients.

Expand

\[
\det\frac{\partial(A,F,G)}{\partial(a,b,c)}+1.
\]

It has 73 distinct source monomials.  The coefficient equations of degree at
most one in the 30 parameters have rank three.  Solving them leaves 27
variables and 64 distinct nonzero equations.

With the free variables ordered first by ambient monomial and then by target
row, the reduced ideal is the unit ideal modulo

\[
32003,\qquad 32009,\qquad 32027.
\]

Singular's modular-over-\(\mathbb Q\) reconstruction

```text
modGB("slimgb",I,1)
```

returns the exact reduced basis

\[
\boxed{\operatorname{GB}(I_{\mathbb Q})=\{1\}}.
\]

Hence the normalized coefficient scheme is empty over \(\mathbb Q\), and
therefore over every characteristic-zero field by the Lefschetz principle
and coefficient descent. \(\square\)

The calculation is reproduced by

```bash
.venv/bin/python scripts/verify_danielewski_one_dicritical_quadratic_obstruction.py
```

## 2. A descent criterion for a displayed Keller map

Let

\[
P(s)=sQ(s),\qquad Q(0)\ne0,
\]

and consider the general Danielewski open-immersion ring

\[
R_Q=k[a,c,s=bc,w=bQ(s)]\subset k[a,b,c].           \tag{4}
\]

Every polynomial in \(k[a,b,c]\) has a unique monomial normal form obtained
by replacing

\[
b^pc^q=
\begin{cases}
b^{p-q}s^q,&p\ge q,\\
c^{q-p}s^p,&q>p.
\end{cases}                                       \tag{5}
\]

Write its positive-\(b\) part as

\[
\sum_{i>0} b^i h_i(a,s).
\]

### Proposition 2.1

A polynomial \(H\in k[a,b,c]\) belongs to \(R_Q\) if and only if

\[
Q(s)^i\mid h_i(a,s)
\qquad\text{for every }i>0.                        \tag{6}
\]

Indeed, the positive-\(b\) normal forms of elements of (4) are exactly

\[
w^i\widetilde h_i(a,s)
=b^iQ(s)^i\widetilde h_i(a,s).
\]

The remaining terms in (5) already lie in \(k[a,c,s]\).  This proves both
directions.

## 3. The foundational counterexample does not descend by permutation

Put \(u=1+xy\) and recall the announced Keller map

\[
\begin{aligned}
K_1&=u^3z+y^2u(4+3xy),\\
K_2&=y+3xu^2z+3xy^2(4+3xy),\\
K_3&=2x-3x^2y-x^3z.
\end{aligned}                                      \tag{7}
\]

For each of the six assignments of \((x,y,z)\) to the Danielewski source
roles \((a,b,c)\), apply (5) to all three coordinates in (7).  For every
positive \(b^i\), take the gcd in \(k[s]\) of the coefficients of
\(h_i(a,s)\) with respect to \(a\).  Condition (6) says that every
irreducible factor of \(Q\) must occur in that gcd with multiplicity at
least \(i\).

In all six assignments the common allowable factor is

\[
\boxed{Q=1}.                                       \tag{8}
\]

Consequently no coordinate permutation places all three coordinates of
(7) in a smooth Danielewski ring (4) with a missing divisor.

This statement is intentionally narrower than a conjugacy theorem.
Nonlinear polynomial source automorphisms, or a different affine
modification chart, are not covered by (8).

The six exact divisibility audits are reproduced by

```bash
.venv/bin/python scripts/audit_foundational_danielewski_descent.py
```

## 4. Consequences

The one-dicritical and two-dicritical targets now have parallel low-degree
frontiers:

\[
\begin{array}{c|c|c}
P(x)&\text{missing planes}&\text{ambient-quadratic contraction}\\ \hline
x(x+1)&1&\text{impossible}\\
x(x^2+1)&2&\text{impossible}.
\end{array}
\]

The first possible one-dicritical contraction has ambient degree at least
three.  Since the relation \(cw=x^2+x\) already has degree two, the
degree-at-most-three quotient has dimension 30.  After translation and
linear-jet normalization, its universal cubic chart has

\[
3(30-1-3)=78
\]

parameters, rather than the 90 parameters of the cubic two-dicritical
chart.

The foundational counterexample shows that affine-linear dependence on one
source coordinate is not by itself an obstruction.  The next useful search
must therefore exploit the Danielewski divisibility and boundary valuation,
not merely suspension degree.  Two concrete continuations remain:

1. construct the 78-variable cubic flux and lowest-coefficient ledger for
   \(x(x+1)\), branching first by boundary pole type;
2. search for filtration-controlled polynomial source automorphisms that
   move the foundational map into some ring \(R_Q\), using Proposition 2.1
   as an exact membership filter.

Either direction directly tests whether the known Keller mechanism can be
realized as a volume-preserving contraction of a smooth Danielewski
completion.
