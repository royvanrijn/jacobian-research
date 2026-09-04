# The septuple-linear Hessian pole--defect ladder

## Status

This note proves `HC4NHM2`: the generic-corank-one Schur packet with
`det(Hess(h5))=ell^7*R2`, where `R2` is squarefree and coprime to `ell`, is
empty. It does **not itself** treat multiplicity eight or nine, a
nonsquarefree quadratic cofactor, lower-Smith boundaries, or unrestricted
`HC(4)`. The first two cases are subsequently closed by `HC4NHM3`.

The nonreduced Hessian--Schur module first turns the possible pole orders
along a septuple line into a six-row finite ladder. Hessian integrability
excludes both degree-one kernel rows and the quadratic-pencil part of the
extremal degree-two row. The resulting conic-kernel branch admits a complete
binary-quintic classification; its boundary kernel and first-normal equations
then exclude every case. The three constant-kernel rows close by direct
determinant faces. Thus all six rows are empty.

Replay the exact polynomial identities with

~~~bash
.venv/bin/python scripts/verify_hc4_direct_septuple_linear_hessian_gate.py
~~~

The module input is
[`HC4_NONREDUCED_HESSIAN_SCHUR_MODULE.md`](HC4_NONREDUCED_HESSIAN_SCHUR_MODULE.md),
and the preceding repeated-line exclusions are in
[`HC4_DIRECT_DOUBLE_LINEAR_HESSIAN_GATE.md`](HC4_DIRECT_DOUBLE_LINEAR_HESSIAN_GATE.md).

Throughout, one may extend the characteristic-zero ground field to its
algebraic closure. Put \(S=K[x,y,z]\), let

\[
C=\operatorname{Hess}(h_5),\qquad d=\nabla s_3,
\tag{0.1}
\]

and assume

\[
\det C=x^7R_2,\qquad x\nmid R_2,\qquad R_2\text{ squarefree},
\tag{0.2}
\]

with \(C\bmod x\) of generic rank two. The Schur face is

\[
\det C\mid d^{\mathsf T}\operatorname{adj}(C)d.
\tag{0.3}
\]

The line \(L=(x=0)\simeq\mathbf P^1\) is assumed essential for the minimal
denominator clearing \(C^{-1}d\).

## 1. The complete pole--defect ladder

At the generic point of \(L\), the local Smith form from `HC4NHM1` shows that
the pole order \(b\) of \(C^{-1}d\) satisfies

\[
1\le b\le\left\lfloor\frac72\right\rfloor=3.
\tag{1.1}
\]

The squarefree residual factor contributes no pole. Hence

\[
P=x^b,\qquad e=x^bC^{-1}d\in S_{b-1}^3.
\tag{1.2}
\]

Let \(\mathcal K_L\subset\mathcal O_L^3\) be the saturated kernel line of
\(C|_L\), and let \(B_L\) be its normalized corank-two defect divisor. The
self-dual determinant formula is

\[
\mathcal K_L^{\otimes2}\simeq\mathcal O_L(B_L-6).
\tag{1.3}
\]

Because a line subbundle of \(\mathcal O_L^3\) has nonpositive degree, write

\[
\mathcal K_L\simeq\mathcal O_L(-\kappa),\qquad \kappa\ge0.
\tag{1.4}
\]

Equations (1.3) and (1.4) give

\[
\boxed{\deg B_L=6-2\kappa.}
\tag{1.5}
\]

The nonzero residue of \(e\) is a section of
\(\mathcal K_L(b-1)=\mathcal O_L(b-1-\kappa)\), so
\(\kappa\le b-1\). This proves the complete ladder:

| pole \(b\) | primitive kernel degree \(\kappa\) | defect length \(\deg B_L\) |
|---:|---:|---:|
| \(1\) | \(0\) | \(6\) |
| \(2\) | \(0\) | \(6\) |
| \(2\) | \(1\) | \(4\) |
| \(3\) | \(0\) | \(6\) |
| \(3\) | \(1\) | \(4\) |
| \(3\) | \(2\) | \(2\) |

If \(\epsilon\) is a primitive generator of the kernel, its entries are
basepoint-free binary forms of degree \(\kappa\), and

\[
e|_L=g\epsilon,\qquad \deg g=b-1-\kappa.
\tag{1.6}
\]

Thus the defect divisor measures exactly how much of the quadratic residue
is common scalar zero and how much is genuine motion of the kernel line.

## 2. The extremal defect-two row

Take the last row:

\[
b=3,\qquad \deg B_L=2,\qquad
\mathcal K_L\simeq\mathcal O_L(-2).
\tag{2.1}
\]

Now \(e|_L\) is a primitive basepoint-free triple of binary quadratics and
defines

\[
\phi_e:\mathbf P^1\longrightarrow\mathbf P^2.
\tag{2.2}
\]

The coordinate span has dimension two or three. Dimension one is impossible
because a single quadratic has a zero over the algebraic closure. Hence
precisely two projective types occur:

1. the coordinates span \(H^0(\mathbf P^1,\mathcal O(2))\), and
   \(\phi_e\) is a nondegenerate Veronese conic; or
2. the coordinates span a basepoint-free pencil, and \(\phi_e\) is a
   degree-two cover of a projective line.

This classifies the kernel morphism. An independent target change in
\(\mathbf P^2\) need not preserve a Hessian presentation, so one cannot
simply substitute \((y^2,yz,z^2)\) or \((y^2,z^2,0)\) into the Hessian
equations.

## 3. Hessian integrability excludes every moving line image

Use the factorial-normalized boundary expansion

\[
h_5=F_5+xG_4+\frac{x^2}{2}H_3
       +\frac{x^3}{6}J_2+\frac{x^4}{24}K_1
       +\frac{x^5}{120}c.
\tag{3.1}
\]

Along \(L\),

\[
C_0:=C|_L=
\begin{pmatrix}
H& G_y&G_z\\
G_y&F_{yy}&F_{yz}\\
G_z&F_{yz}&F_{zz}
\end{pmatrix}.
\tag{3.2}
\]

Suppose a primitive kernel generator has positive degree
\(\kappa\in\{1,2\}\) and its entries span a plane \(W\subset K^3\). Put
\(T_L=\langle\partial_y,\partial_z\rangle\). For \(\kappa=1\), a
basepoint-free triple of linear forms necessarily has this form. For
\(\kappa=2\), this is the pencil alternative from Section 2.

### 3.1 A transverse kernel plane drops the boundary rank

If \(W\ne T_L\), a coordinate change preserving \(L\) puts

\[
W=\langle\partial_x,\partial_y\rangle,\qquad e|_L=(a,b,0),
\tag{3.3}
\]

where \(a,b\) are independent basepoint-free forms of degree \(\kappa\).
The last two rows
of \(C_0e=0\) give, where \(a\ne0\),

\[
\nabla G=-R\,\partial_y\nabla F,\qquad R=\frac ba.
\tag{3.4}
\]

Curl-freeness gives

\[
R_zF_{yy}=R_yF_{yz}.
\tag{3.5}
\]

The ratio \(R\) is nonconstant. Hence, generically,
\((F_{yy},F_{yz})\) and \((R_y,R_z)\) are proportional. Euler's identities

\[
yR_y+zR_z=0,\qquad yF_{yy}+zF_{yz}=4F_y
\tag{3.6}
\]

force \(F_y=0\). Equations (3.4) and the first row of \(C_0e=0\) then give
\(G=H=0\). Thus \(C_0\) has rank at most one, a contradiction.

### 3.2 A tangent kernel plane gives only a simple line

If \(W=T_L\), write

\[
e|_L=(0,b,c)
\tag{3.7}
\]

with \(b,c\) independent and basepoint-free. The tangent rows give

\[
\operatorname{Hess}(F)(b,c)^{\mathsf T}=0.
\tag{3.8}
\]

The binary zero-Hessian theorem makes a nonzero \(F\) a fifth power. Its
Hessian has a constant kernel direction, incompatible with the moving
pencil \((b,c)\). Hence \(F=0\). The remaining equation is

\[
G_yb+G_zc=0.
\tag{3.9}
\]

Let \(A=\gcd(G_y,G_z)\). The primitive syzygy is
\((-G_z/A,G_y/A)\), of degree \(3-\deg A\). Hence

\[
\deg A=3-\kappa.
\tag{3.9a}
\]

For a binary quartic, \(\deg A\) is its total repeated-root excess. Thus
\(\kappa=1\) gives the root partitions \(3+1\) and \(2+2\), while
\(\kappa=2\) gives \(2+1+1\). In every case \(G\) has at least two distinct
roots.

The coefficient of \(x\) in the full Hessian determinant is independent of
all lower binary jets. Jacobi's first-variation formula gives

\[
\begin{aligned}
[x]\det C
&=-G_z^2G_{yy}+2G_yG_zG_{yz}-G_y^2G_{zz}\\
&=-\frac43G\det\operatorname{Hess}(G).
\end{aligned}
\tag{3.10}
\]

The last equality follows from Euler's identity. A quartic of root type
\(2+1+1\) has nonzero Hessian, so (3.10) is nonzero. Therefore the
determinant line has multiplicity exactly one.

Both positions of the kernel plane are impossible when
\(x^2\mid\det C\). Hence

\[
\boxed{
\text{A positive-degree kernel whose image is a projective line is
incompatible with a repeated Hessian line.}
}
\tag{3.11}
\]

For \(\kappa=1\), the image is always a line. Consequently both
defect-four ladder rows

\[
(b,\kappa,\deg B_L)=(2,1,4),\qquad(3,1,4)
\tag{3.12}
\]

are empty. For \(\kappa=2\), only the nondegenerate conic alternative
survives this section; Sections 4--6 exclude it as well.

## 4. The conic boundary and first-normal gate

Write the conic kernel in the original Hessian coordinates as

\[
e|_L=(a,b,c),\qquad a,b,c\in K[y,z]_2.
\tag{4.1}
\]

The entries are linearly independent, so \(a\ne0\). Since \(C_0\) has rank
two and defect degree two, a nonzero quadratic \(\lambda\) satisfies

\[
\operatorname{adj}(C_0)=\lambda ee^{\mathsf T}.
\tag{4.2}
\]

The \(xx\) cofactor gives the square-Hessian equation

\[
\boxed{\det\operatorname{Hess}(F)=\lambda a^2.}
\tag{4.3}
\]

The first normal derivative of \(C\) is

\[
C_1=
\begin{pmatrix}
J&H_y&H_z\\
H_y&G_{yy}&G_{yz}\\
H_z&G_{yz}&G_{zz}
\end{pmatrix}.
\tag{4.4}
\]

Using (4.2),

\[
[x]\det C
=\lambda\left(
a^2J+2a(bH_y+cH_z)
 +(b,c)\operatorname{Hess}(G)(b,c)^{\mathsf T}
\right).
\tag{4.5}
\]

A repeated line therefore requires

\[
\boxed{
a^2\mid
2a(bH_y+cH_z)
 +(b,c)\operatorname{Hess}(G)(b,c)^{\mathsf T}.
}
\tag{4.6}
\]

When (4.6) holds, (4.5) determines \(J\) uniquely. The septuple condition
then consists of \([x^2]\det C,\ldots,[x^6]\det C=0\), followed by a
nonzero degree-two coefficient at \(x^7\).

Only afterward should one impose the cleared Schur equations

\[
Ce=x^3d,\qquad d=\nabla s_3,\qquad d^{\mathsf T}e=x^3\alpha_1.
\tag{4.7}
\]

The last equation is equivalent to the original Schur divisibility because
\(\operatorname{adj}(C)d=x^4R_2e\).

## 5. Repeated-root sieve for the binary quintic

Equation (4.3) can already be classified on every repeated-root stratum of
\(F\). In an affine coordinate \(t\), write \(F=z^5f(t)\). Up to a nonzero
scalar, its binary Hessian is

\[
5ff''-4(f')^2.
\tag{5.1}
\]

If \(f=t^m u(t)\), with \(u(0)\ne0\) and \(1\le m\le4\), the initial term
of (5.1) is

\[
m(m-5)u(0)^2t^{2m-2}.
\tag{5.2}
\]

Thus a root of multiplicity \(m\) in \(F\) has Hessian multiplicity exactly
\(2m-2\). It follows immediately that the root partitions

\[
4+1,\qquad 3+2,\qquad 3+1+1,\qquad 2+2+1
\tag{5.3}
\]

have a square Hessian divisor of degree at least four. The pure fifth power
has zero Hessian and is incompatible with the nonzero \(\lambda\) in (4.3).

The only repeated-root partition not decided by (5.2) is
\(2+1+1+1\). Normalize it to

\[
F_t=u^2v(u-v)(u-tv),\qquad t\ne0,1.
\tag{5.4}
\]

Its Hessian is \(u^2Q_t(u,v)\), where \(Q_t\) is quartic. A second linear
square occurs exactly when \(Q_t\) has a repeated root. Exact elimination
gives

\[
\operatorname{disc}(Q_t)
\doteq t^2(t-1)^2P(t),
\tag{5.5}
\]

where

\[
P(t)=4t^6-12t^5-t^4+22t^3-t^2-12t+4.
\tag{5.6}
\]

The sextic is the anharmonic-invariant equation

\[
P(t)=4(t^2-t+1)^3-25t^2(1-t)^2.
\tag{5.7}
\]

Consequently all six roots form one orbit under the transformations
generated by \(t\mapsto1-t\) and \(t\mapsto1/t\). Equivalently, the unique
exceptional orbit is

\[
\boxed{
\frac{(t^2-t+1)^3}{t^2(1-t)^2}=\frac{25}{4}.
}
\tag{5.8}
\]

Hence the conic square-Hessian gate has the following repeated-root input:
the four partitions in (5.3), plus one exceptional
\(2+1+1+1\) orbit. This is only a necessary boundary sieve: one must still
solve \(C_0e=0\) with the actual normal component \(a\) and apply (4.6).

### 5.1 The squarefree stratum is one orbit

The squarefree calculation is also finite. Normalize five distinct roots to

\[
F_{s,t}=uv(u-v)(u-sv)(u-tv),
\qquad
D=s\,t(s-1)(t-1)(s-t)\ne0.
\tag{5.9}
\]

On the affine chart \(v=1\), write a possible quadratic square factor and
residual quadratic as

\[
q=u^2+Au+B,\qquad
\lambda=-16u^2+L_1u+L_0.
\tag{5.10}
\]

This chart loses no projective orbit: among five distinct roots of \(F\),
one can choose the root sent to infinity away from the two roots of \(q\).
Comparing the seven coefficients in

\[
\det\operatorname{Hess}(F_{s,t})=q^2\lambda
\tag{5.11}
\]

first eliminates \(L_1,L_0\). Saturating the four remaining equations by
\(D\) and eliminating \(A,B\) gives the exact ideal

\[
\begin{aligned}
E_1={}&s^3-s^2t-st^2+t^3-s^2+2st-t^2-s-t+1,\\
E_2={}&s^2t^2+2st^3-t^4-s^2t-4st^2+t^3
        +s^2-4st+2t^2+2s+t-1,\\
E_3={}&3st^4-t^5-6st^3+t^4+4s^2t-7st^2+3t^3
        -2s^2+6st+t^2-s-3t+1,\\
E_4={}&t^6-3t^5-2t^4+9t^3-2t^2-3t+1.
\end{aligned}
\tag{5.12}
\]

Its reduced lexicographic basis consists of \(E_4\) and

\[
\begin{aligned}
S={}&2s^2-2st^5+5st^4+8st^3-17st^2-8st+5s\\
   &\quad+t^5-3t^4-3t^3+9t^2+3t-3.
\end{aligned}
\tag{5.13}
\]

The terminal polynomial factors as

\[
E_4=(t^2-3t+1)(t^2-t-1)(t^2+t-1).
\tag{5.14}
\]

Modulo these three factors, \(S/2\) is respectively

\[
s^2-t,\qquad s^2-2st+t,\qquad s^2-2s+t.
\tag{5.15}
\]

Thus the saturated squarefree scheme is reduced of degree twelve. The six
possible values of \(t\) form one anharmonic orbit, characterized by

\[
\frac{(t^2-t+1)^3}{t^2(1-t)^2}=8.
\tag{5.16}
\]

To identify the projective orbit, take
\(\varphi^2-\varphi-1=0\). One representative is

\[
\{0,\infty,1,\varphi,\varphi^2\}.
\tag{5.17}
\]

The two roots \(s=\pm\varphi\) in the same \(t=\varphi^2\) chart are
projectively equivalent: \(x\mapsto x/(x-\varphi)\) sends the first
configuration to the second after renormalizing the three marked roots.
The Möbius transformation

\[
T(x)=\frac{1}{1+(\varphi-2)x}
\tag{5.18}
\]

cycles

\[
\infty\longmapsto0\longmapsto1\longmapsto\varphi
\longmapsto\varphi^2\longmapsto\infty.
\tag{5.19}
\]

Over the algebraic closure, an order-five element of
\(\operatorname{PGL}_2\) is diagonalizable. Conjugating \(T\) to
\(x\mapsto\zeta_5x\) makes (5.17) one nonfixed orbit, so its binary quintic
is projectively equivalent to \(u^5+v^5\). Conversely,

\[
\det\operatorname{Hess}(u^5+v^5)
=400u^3v^3=(uv)^2(400uv).
\tag{5.20}
\]

Therefore

\[
\boxed{
\text{The only squarefree binary quintic with a quadratic square divisor
in its Hessian is the orbit of }u^5+v^5.
}
\tag{5.21}
\]

This orbit is nevertheless incompatible with the primitive conic kernel.
For \(F=u^5+v^5\), equation (4.3) reads

\[
400u^3v^3=\lambda a^2.
\tag{5.22}
\]

Because \(a\) is quadratic, unique factorization forces
\(a\doteq uv\). Write \(e=(uv,b,c)\). The last two rows of \(C_0e=0\)
become

\[
20u^3b=-uvG_u,\qquad
20v^3c=-uvG_v.
\tag{5.23}
\]

Polynomiality gives

\[
u^2\mid G_u,\qquad v^2\mid G_v.
\tag{5.24}
\]

For a binary quartic this forces

\[
G=\alpha u^4+\beta v^4.
\tag{5.25}
\]

Equation (5.23) then gives

\[
e=uv\left(1,-\frac{\alpha}{5},-\frac{\beta}{5}\right).
\tag{5.26}
\]

All three entries have the common factor \(uv\), contradicting the
primitive basepoint-free kernel required on the defect-two row. Hence

\[
\boxed{
\text{The squarefree binary-quintic boundary of the extremal conic row is
empty.}
}
\tag{5.27}
\]

Together with (5.3) and (5.8), the surviving binary input is therefore
entirely repeated-root. The first-normal gate and higher determinant
coefficients remain to be imposed on those repeated-root strata.

## 6. Completion of the extremal defect-two row

It remains to combine the binary list with the full boundary kernel equation.
The cases with at most three distinct roots have no moduli. Direct solution
of

\[
\operatorname{Hess}(F)(b,c)^{\mathsf T}
=-a\nabla G
\tag{6.1}
\]

gives the following table. The displayed families are the complete linear
solution spaces.

| root type | \(F\) | possible \(a\) | tangent solution |
|---|---|---|---|
| \(3+1+1\) | \(u^3v(u-v)\) | \(u^2\) | \((b,c)=(p u^2,q u^2)\) |
| \(2+2+1\) | \(u^2v^2(u-v)\) | \(uv\) | \((b,c)=(p uv,q uv)\) |
| \(3+2\) | \(u^2v^3\) | \(uv\) | \((b,c)=(p uv,q uv)\) |

In all three rows, \(a,b,c\) have a common quadratic factor and therefore
do not define a primitive conic kernel.

For \(F=u^2v^3\), there is a second square divisor \(a=v^2\). Its complete
tangent solution is

\[
\begin{aligned}
G&=-\frac{uv}{2}(4pv^2+3qu^2+6ruv),\\
b&=\frac{4pv^2-3qu^2}{4},\\
c&=v(qu+rv).
\end{aligned}
\tag{6.2}
\]

The first row of \(C_0e=0\) requires
\(a\mid G_ub+G_vc\). The remainder modulo \(v^2\) is

\[
\frac{15}{8}q^2u^4v.
\tag{6.3}
\]

Thus \(q=0\), after which \(a,b,c\) are all divisible by \(v^2\).
This closes the \(3+2\) partition.

### 6.1 The \(4+1\) partition reaches the first-normal gate

Take \(F=uv^4\). Its Hessian determinant is \(-16v^6\), so
\(a=v^2\). The complete solution of \(C_0e=0\) is

\[
\begin{aligned}
G={}&-v^2(p_0v^2+2p_1u^2+4p_2uv),\\
b={}&p_0v^2-2p_1u^2,\\
c={}&v(p_1u+p_2v),\\
H={}&4(p_1u+p_2v)(2p_0v^2-p_1u^2+3p_2uv).
\end{aligned}
\tag{6.4}
\]

The kernel is primitive and spans all binary quadratics exactly when
\(p_1\ne0\). The first-normal numerator from (4.6) has remainder

\[
60p_1^3u^4v^2\pmod{v^4}.
\tag{6.5}
\]

Divisibility by \(a^2=v^4\) forces \(p_1=0\), precisely the
nonprimitive/nonconic locus. Hence the \(4+1\) partition is empty.

The earlier calibration is the point
\((p_0,p_1,p_2)=(0,-1/2,0)\). Namely,

The conic boundary is nonempty before the first-normal condition. For binary
variables \(u,v\), take

\[
F=uv^4,\qquad G=u^2v^2,\qquad H=-u^3.
\tag{6.6}
\]

Then

\[
C_0=
\begin{pmatrix}
-u^3&2uv^2&2u^2v\\
2uv^2&0&4v^3\\
2u^2v&4v^3&12uv^2
\end{pmatrix}
\tag{6.7}
\]

has

\[
e=(-2v^2,-2u^2,uv)^{\mathsf T},\qquad
\operatorname{adj}(C_0)=-4v^2ee^{\mathsf T}.
\tag{6.8}
\]

The entries of \(e\) span all binary quadratics. However, the numerator in
(4.6) is

\[
2a(bH_u+cH_v)
 +(b,c)\operatorname{Hess}(G)(b,c)^{\mathsf T}
=-30u^4v^2,
\tag{6.9}
\]

which is not divisible by \(a^2=4v^4\). This point does not extend even to a
double determinant line. It calibrates the conic stratum and the nonvacuity
of the first-normal gate without claiming a septuple example.

### 6.2 The exceptional \(2+1+1+1\) orbit

The orbit (5.8) has the rational normal form

\[
F=v^2(5u^3+30uv^2+8v^3).
\tag{6.10}
\]

Its binary Hessian factors as

\[
\det\operatorname{Hess}(F)
=-600v^2(u-2v)^2(u^2+4uv+6v^2).
\tag{6.11}
\]

Thus the unique quadratic square divisor is

\[
a=v(u-2v).
\tag{6.12}
\]

The complete tangent solution is

\[
\begin{aligned}
G={}&\frac{5v}{2}
 \left(3pu^2v+6pv^3+2qu^3+24quv^2+8qv^3\right),\\
b={}&-\frac p2v(u-2v),\\
c={}&-\frac q2v(u-2v).
\end{aligned}
\tag{6.13}
\]

Again \(a,b,c\) share the quadratic factor \(v(u-2v)\), so the kernel is
not primitive. Section 5 already excluded the squarefree Fermat orbit.
Every binary quintic from the complete square-Hessian classification has
now been removed.

> **Extremal septuple-line exclusion.** Under (0.1)--(0.3), the
> pole-three/defect-two row \((b,\kappa,\deg B_L)=(3,2,2)\) is empty.

Together with (3.12), every positive-degree kernel row is empty. The
septuple-line ladder reduces to

\[
\boxed{
(b,\kappa,\deg B_L)=(1,0,6),\ (2,0,6),\ (3,0,6).
}
\tag{6.14}
\]

These are precisely the constant-kernel rows.

## 7. The constant-kernel rows

Let \(\epsilon\in K^3\) be the constant primitive kernel of \(C_0\).
There are two orbits under linear coordinate changes preserving \(L\).

### 7.1 A transverse constant kernel

If \(\epsilon\notin T_L\), normalize \(\epsilon=\partial_x\). Equation
\(C_0\epsilon=0\) gives

\[
G=H=0.
\tag{7.1}
\]

Thus

\[
h_5=F+\frac{x^3}{6}J+\frac{x^4}{24}K+\frac{x^5}{120}c,
\tag{7.2}
\]

and generic rank two along \(L\) says
\(\det\operatorname{Hess}(F)\ne0\). Successive determinant coefficients are

\[
\begin{aligned}
[x]\det C&=J\det\operatorname{Hess}(F),\\
[x^2]\det C\big|_{J=0}
 &=\frac K2\det\operatorname{Hess}(F),\\
[x^3]\det C\big|_{J=K=0}
 &=\frac c6\det\operatorname{Hess}(F).
\end{aligned}
\tag{7.3}
\]

Divisibility by \(x^7\) forces \(J=K=c=0\), after which
\(\det C=0\). Hence the transverse constant-kernel branch is empty.

### 7.2 A tangent constant kernel

If \(\epsilon\in T_L\), normalize \(\epsilon=\partial_z\). The boundary
kernel equations give

\[
F=\alpha y^5,\qquad G=\beta y^4,
\tag{7.4}
\]

while

\[
H=h_0y^3+h_1y^2z+h_2yz^2+h_3z^3
\tag{7.5}
\]

is initially arbitrary. Write

\[
J=j_0y^2+j_1yz+j_2z^2,\qquad
K=k_0y+k_1z.
\tag{7.6}
\]

The induced rank-two determinant is

\[
20\alpha y^3H-16\beta^2y^6.
\tag{7.7}
\]

There are two cases.

If \(\alpha\ne0\), the \(x^2\)-coefficient successively gives

\[
h_3=h_2=h_1=0.
\tag{7.8}
\]

Generic rank two is then

\[
5\alpha h_0-4\beta^2\ne0.
\tag{7.9}
\]

The next faces are

\[
\begin{aligned}
[x^3]\det C&=\frac43j_2y^6(5\alpha h_0-4\beta^2),\\
[x^4]\det C\big|_{j_2=0}&=-5\alpha j_1^2y^5,\\
[x^6]\det C\big|_{j_2=j_1=0}&=-\frac59\alpha k_1^2y^3.
\end{aligned}
\tag{7.10}
\]

Thus \(j_2=j_1=k_1=0\), and then

\[
[x^7]\det C=0.
\tag{7.11}
\]

So the determinant cannot have exact \(x\)-multiplicity seven.

If \(\alpha=0\), generic rank two forces \(\beta\ne0\). The \(x^2\) and
\(x^3\) faces give

\[
h_2=h_3=0,\qquad
j_2=-\frac{3h_1^2}{4\beta}.
\tag{7.12}
\]

The coefficient of \(x^4y^4z\) is \(15h_1^3\), so \(h_1=0\) and
\(j_2=0\). The next face is

\[
[x^5]\det C=-\frac73\beta j_1^2y^4,
\tag{7.13}
\]

hence \(j_1=0\). Finally,

\[
\boxed{
[x^7]\det C=-\frac13\beta k_1^2y^2.
}
\tag{7.14}
\]

If \(k_1=0\), the multiplicity is greater than seven. If \(k_1\ne0\), the
residual quadratic is a square, contrary to (0.2). Thus the tangent
constant-kernel branch is empty as well.

Combining Sections 3, 6, and 7 proves:

> **Theorem `HC4NHM2` -- Exact septuple-linear Hessian--Schur exclusion.**
> Under (0.1)--(0.3), with \(C\bmod x\) of generic rank two and \(x\)
> essential for the minimal Schur denominator, no packet with
> \(\det C=x^7R_2\), \(R_2\) squarefree and \(x\nmid R_2\), exists.

The proof actually eliminates the Hessian boundary before the curl and scalar
Schur equations are needed. It does not treat linear multiplicity eight or
nine, a nonsquarefree residual quadratic, or a lower-Smith boundary.

## 8. Historical next calculation and current boundary

The next repeated-line targets from this theorem are

\[
\det C=x^8\ell,\qquad \det C=x^9,
\tag{8.1}
\]

on the generic-corank-one stratum, together with the separate generic
corank-two/three Smith strata. The generic-corank-one targets in (8.1) are
closed in
[`HC4_OCTUPLE_NONUPLE_LINEAR_HESSIAN_GATE.md`](HC4_OCTUPLE_NONUPLE_LINEAR_HESSIAN_GATE.md).
The exact-septuple squarefree-cofactor packet is closed here. The septuple
line with nonsquarefree quadratic cofactor is the `x^7*y^2` incidence later
classified by `HC4NHM4` and excluded from four-variable prolongation by
`HC4NHM5`. Lower-Smith branches remain outside this theorem.
