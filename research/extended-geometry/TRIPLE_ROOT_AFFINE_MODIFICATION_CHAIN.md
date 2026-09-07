# The unified triple-root affine-modification chain

The square-oriented and primitive-multiplicity targets belong to one
residue-crepant affine-modification chain which resolves to the affine
source.  Singularities and finite-field defects change only by moving,
filling, or merging boundary components.  The final affine-three target is
the source itself and carries no dicritical divisor.

## 1. Five stages

Put

\[
 f(x,y)=x\,y(x+1)(x-y+1).
\]

The square-oriented stage is

\[
 T_{-1}=\{S^2=Cf(x,y)\}
 \subset\mathbb A^4_{x,y,C,S}.                       \tag{1}
\]

Adjoining `u=S/C` and writing `v=C` gives the strict transform

\[
 T_0=\{u^2v=f(x,y)\}.                                \tag{2}
\]

Indeed, the morphism

\[
T_0\longrightarrow T_{-1},\qquad C=v,\quad S=uv    \tag{3}
\]

satisfies `(uv)^2=v(u^2v)=vf`.

Now adjoin successively

\[
 a=\frac yu,\qquad c=\frac ua,\qquad b=\frac xc.    \tag{4}
\]

The remaining strict transforms are

\[
\begin{aligned}
T_1&=
\{uv=a\,x(x+1)(x-au+1)\}
\subset\mathbb A^4_{x,a,u,v},\\
T_2&=
\{cv=x(x+1)(x-a^2c+1)\}
\subset\mathbb A^4_{x,a,c,v},\\
T_3&=\mathbb A^3_{a,b,c}.
\end{aligned}                                       \tag{5}
\]

Their morphisms are

\[
T_1\longrightarrow T_0,\qquad y=au,                 \tag{6}
\]

\[
T_2\longrightarrow T_1,\qquad u=ac,                 \tag{7}
\]

\[
T_3\longrightarrow T_2,\qquad
x=bc,\quad
v=b(bc+1)(bc-a^2c+1).                               \tag{8}
\]

The composite `T_3 -> T_0` is the
[primitive-multiplicity map](PRIMITIVE_MULTIPLICITY_TRIPLE_ROOT_MAP.md);
the composite `T_3 -> T_(-1)` is the simplified form of the
[square-oriented map](AFFINE_SOURCE_TRIPLE_ROOT_COX_MAP.md).

## 2. Residue forms telescope

Use

\[
\Omega_{-1}
=\frac{dx\wedge dC\wedge dy}{2S}.
\]

Differentiating `u^2v=f(x,y)` and wedging with `dx,dy` gives

\[
(3)^*\Omega_{-1}
=\frac{dx\wedge dy\wedge du}{u^2}
=:\Omega_0.                                         \tag{9}
\]

The remaining pullbacks are

\[
(6)^*\Omega_0
=\frac{dx\wedge da\wedge du}{u}
=:\Omega_1,                                         \tag{10}
\]

\[
(7)^*\Omega_1
=\frac{dx\wedge da\wedge dc}{c}
=:\Omega_2,                                         \tag{11}
\]

\[
(8)^*\Omega_2
=-da\wedge db\wedge dc.                             \tag{12}
\]

Thus the constant Jacobian is the telescoping product of four primitive
affine modifications.

## 3. Singular loci

The square target has singular locus

\[
\{C=S=f=0\}
\ \cup\
\{S=0,\ f=\partial_xf=\partial_yf=0\}.               \tag{13}
\]

It contains four branch lines at `C=0` and three `C`-lines over the
arrangement intersection points `(0,0),(0,1),(-1,0)`.

The later stages have:

\[
\operatorname{Sing}(T_0)
=
\bigcup_{p\in\{(0,0),(0,1),(-1,0)\}}
\{(x,y,u)=p_1,p_2,0\},                              \tag{14}
\]

with `v` free;

\[
\operatorname{Sing}(T_1)
=
\{x=a=u=v=0\}
\cup
\{x=-1,\ u=v=0\},                                  \tag{15}
\]

where `a` is free in the second component;

\[
\operatorname{Sing}(T_2)
=\{x=-1,\ c=v=0\},                                  \tag{16}
\]

where `a` is free.  Finally, `T_3` is smooth.

## 4. Dicritical tradeoff

The numbers of distinct affine-target dicritical components are

\[
\boxed{4,\quad3,\quad1,\quad1,\quad0}               \tag{17}
\]

for `T_(-1),T_0,T_1,T_2,T_3`.

At `T_(-1)`, the components `C=0`, `x=0`, `x+1=0`, and
`x-y+1=0` are dicritical, while `y=0` is finite.  Adjoining `u=S/C`
moves the `C=0` component to infinity of `T_0`.

At `T_0`, the components `x=0`, `x+1=0`, and `x-y+1=0` are dicritical.
Adjoining `a=y/u` fills `x=0` and makes `x+1` and `x-y+1` coincide along
`u=0`.  Hence only the doubled `x+1=0` component remains dicritical on
`T_1`.

That component remains on `T_2`.  The last coordinate `b=x/c` fills it,
leaving the identity target `T_3=A^3`.

## 5. Exact finite-field behavior

For every finite field,

\[
\begin{array}{c|c|c}
\text{stage}&\#T_i(\mathbb F_q)&
\#T_i(\mathbb F_q)-q^3\\ \hline
T_{-1}&q^3&0,\\
T_0&q(q-1)(q+4)&q(3q-4),\\
T_1&q(q^2+2q-2)&2q(q-1),\\
T_2&q^2(q+1)&q^2,\\
T_3&q^3&0.
\end{array}                                         \tag{18}
\]

The square target already has the correct total point count, but its map
has large missing and multiple fibers.  The primitive modification
simplifies the boundary at the cost of positive point-count excess.  That
excess then decreases and vanishes only at the final identity stage.

### Theorem 5.1

Within the canonical chain (4), no stage simultaneously has an affine-three
target and a nonzero dicritical boundary.  Correct point count alone is
insufficient: the two `q^3` stages are respectively the nonbijective square
target and the dicritical-free identity.

Any successful alternative must leave this chain before adjoining `b=x/c`
and cancel the remaining boundary defect through new geometry rather than
by filling the last dicritical component.

There is one optimal balanced last step.  Adjoining
`w=x(x+1)/c` gives the
[smooth Danielewski target](SMOOTH_DANIELEWSKI_TRIPLE_ROOT_MAP.md).
It preserves polynomiality and one dicritical component, but its class is
`L^3+L^2`, not `L^3`.

## 6. Reproduction

Run

```bash
.venv/bin/python scripts/verify_triple_root_affine_modification_chain.py
```

The checker verifies all five stages, four ring maps, telescoping residue
forms, singular loci, dicritical counts, and finite-field formulas.
