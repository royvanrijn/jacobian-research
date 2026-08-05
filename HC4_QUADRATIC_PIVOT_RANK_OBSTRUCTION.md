# Quadratic-pivot rank obstruction for HC4

## Status

This note begins the nonlinear scalar-pivot branch left open by HC4RSD7.
It completely removes quadratic pivots whose constant Hessian has rank at
least three from the singular-pencil reverse-Schur programme.

> **Theorem HC4RSD8 (quadratic-pivot rank obstruction).** Let \(K\) have
> characteristic zero, let \(x=(x_1,\ldots,x_4)\), and put
>
> \[
> \Phi(t,x)=\frac{\lambda}{2}t^2+tA(x)+B(x),
> \qquad \deg A\leq 2.
> \tag{0.1}
> \]
>
> Suppose
>
> \[
> \det\operatorname{Hess}_{t,x}\Phi=c\in K^\times,
> \qquad
> \det\operatorname{Hess}_x(B+sA)=0.
> \tag{0.2}
> \]
>
> Then
>
> \[
> \operatorname{rank}\operatorname{Hess}A\leq2,
> \qquad
> dA|_{\ker\operatorname{Hess}A}\ne0.
> \tag{0.3}
> \]
>
> Consequently, if \(r=\operatorname{rank}\operatorname{Hess}A\), constant
> affine coordinates put the pivot in the form
>
> \[
> A=w+\frac12u^{\mathsf T}Q_ru,
> \qquad r\leq2,
> \quad \det Q_r\ne0.
> \tag{0.4}
> \]
>
> In rank two, with passive coordinates \((z,w)\),
>
> \[
> \det\operatorname{Hess}_{z,w}B=0.
> \tag{0.5}
> \]
>
> Its rank-zero stratum has a parameter-independent primitive pencil-kernel
> line in the fixed passive plane and is already covered by HC4RSD5.
> Therefore every genuinely new rank-two component has passive Hessian rank
> exactly one.
>
> The bound is sharp: rank-two examples satisfying (0.2) exist. The theorem
> does not classify their moving kernel lines.

The exact block identities are replayed by
[scripts/verify_hc4_quadratic_pivot_rank_obstruction.py](scripts/verify_hc4_quadratic_pivot_rank_obstruction.py),
which writes
[artifacts/generated-results/hc4_quadratic_pivot_rank_obstruction.json](artifacts/generated-results/hc4_quadratic_pivot_rank_obstruction.json).

## 1. Rank four is immediate

Put

\[
Q=\operatorname{Hess}A,
\qquad
M(s,x)=\operatorname{Hess}B+sQ.
\tag{1.1}
\]

The coefficient of \(s^4\) in \(\det M\) is \(\det Q\). Thus the singular
pencil in (0.2) immediately excludes \(\operatorname{rank}Q=4\).

## 2. The rank-three block

Assume \(\operatorname{rank}Q=3\). After a constant linear change of the
four \(x\)-coordinates, write \(x=(u,z)\), with
\(u=(u_1,u_2,u_3)\), so

\[
Q=\begin{pmatrix}Q_3&0\\0&0\end{pmatrix},
\qquad \det Q_3\ne0.
\tag{2.1}
\]

Write the pencil Hessian as

\[
M=
\begin{pmatrix}
K&d\\ d^{\mathsf T}&e
\end{pmatrix},
\qquad
K=sQ_3+H,
\tag{2.2}
\]

where \(H,d,e\) are polynomial in \(x\). The coefficient of \(s^3\) in
\(\det M\) is

\[
e\det Q_3.
\tag{2.3}
\]

The singular-pencil equation therefore forces \(e=0\). Equivalently,
\(D_z^2B=0\), so \(B\) is affine in the null direction \(z\). With \(e=0\),

\[
\det M=-d^{\mathsf T}\operatorname{adj}(K)d=0.
\tag{2.4}
\]

## 3. The odd-degree square contradiction

Split the affine gradient of the quadratic pivot as

\[
\nabla A=(g,\gamma),
\qquad g\in K[x]^3,
\quad \gamma\in K.
\tag{3.1}
\]

The full parent Hessian is

\[
\mathcal H=
\begin{pmatrix}
\lambda&g^{\mathsf T}&\gamma\\
g&K&d\\
\gamma&d^{\mathsf T}&0
\end{pmatrix}.
\tag{3.2}
\]

Define

\[
\begin{aligned}
E&=d^{\mathsf T}\operatorname{adj}(K)d,\\
P&=\gamma\det K-g^{\mathsf T}\operatorname{adj}(K)d,\\
G&=\lambda\det K-g^{\mathsf T}\operatorname{adj}(K)g.
\end{aligned}
\tag{3.3}
\]

Exact block elimination, cleared of the denominator \(\det K\), gives

\[
\boxed{
(\det K)\det\mathcal H-G\det M+P^2=0.
}
\tag{3.4}
\]

Using (0.2) and \(\det M=0\), this becomes

\[
P^2=-c\det K.
\tag{3.5}
\]

But \(\det K=\det(sQ_3+H)\) has degree exactly three in \(s\), with leading
coefficient \(\det Q_3\ne0\). The nonzero polynomial on the left of (3.5)
has even \(s\)-degree, while the right side has degree three. This is
impossible. Hence rank three is excluded, proving the rank bound in (0.3).

Notice that the argument uses neither a constant kernel direction nor the
Piola classification. It therefore closes rank three even when the kernel
line of \(M(s,x)\) moves with \(x\) and \(s\).

## 4. The surviving-pivot normal form

The singular pencil and the nonzero parent determinant force the generic
rank of \(M\) to be three. Moreover,

\[
c=-\nabla A^{\mathsf T}\operatorname{adj}(M)\nabla A.
\tag{4.1}
\]

Thus the entries of \(\nabla A\) generate the unit ideal: the right side of
(4.1) is already an explicit polynomial combination of them equal to a
nonzero constant.

Write \(A=\tfrac12x^{\mathsf T}Qx+a^{\mathsf T}x+a_0\). If \(a\) vanished
on \(\ker Q\), then \(a\in\operatorname{im}Q\), and the affine system
\(Qx+a=0\) would have a solution. Every entry of \(\nabla A\) would vanish
there, contradicting the unit ideal. Hence \(a|_{\ker Q}\ne0\).

Choose \(w\in\ker Q\) with \(D_wA=1\), remove the active linear part by a
translation, and choose the remaining null coordinates inside the kernel of
\(a|_{\ker Q}\). This gives (0.4). For \(r=0\), the pivot is affine and is
already closed by HC4RSD7. The genuinely nonlinear quadratic frontier is
therefore exactly \(r=1,2\) in this normal form.

## 5. Sharp rank-two calibration

Let \(x=(x,y,z,w)\) and take

\[
A=yz+w,
\qquad
B=xz+\frac12y^2.
\tag{5.1}
\]

Then \(\operatorname{rank}\operatorname{Hess}A=2\), while

\[
\det\operatorname{Hess}(B+sA)=0,
\qquad
\operatorname{rank}_{K(s)}\operatorname{Hess}(B+sA)=3,
\tag{5.2}
\]

and for every \(\lambda\),

\[
\det\operatorname{Hess}
\left(\frac{\lambda}{2}t^2+tA+B\right)=1.
\tag{5.3}
\]

This example has the fixed kernel direction \(\partial_w\), so it belongs
to the already closed cone-pencil stratum HC4RSD1. It shows why the rank
argument itself cannot be pushed below two.

## 6. Revised quadratic frontier

The singular-pencil quadratic-pivot branch is now reduced to:

1. \(\operatorname{rank}\operatorname{Hess}A=1\) or \(2\);
2. removal of the common-constant-kernel cases already covered by
   HC4RSD1--HC4RSD5;
3. classification of the residual \(x\)-moving corank-one kernel lines;
4. collision and compactified-gradient tests only on those residual
   components.

For rank two, use the normal form

\[
A=w+\frac12u^{\mathsf T}Q_2u
\]

with passive coordinates \((z,w)\). The coefficient of \(s^2\) in the
singular pencil determinant is

\[
(\det Q_2)\det\operatorname{Hess}_{z,w}B.
\]

Thus the passive binary Hessian of \(B\) is singular over \(K(u)\). Its
rank-zero part is already closed. Indeed, if
\(\operatorname{Hess}_{z,w}B=0\), then

\[
B=z\,b(u)+w\,c(u)+D_0(u).
\]

Let \(J=(\nabla b,\nabla c)\). In active/passive block coordinates the
pencil has form

\[
M=\begin{pmatrix}K&J\\J^{\mathsf T}&0\end{pmatrix},
\qquad
\det M=(\det J)^2.
\]

Pencil singularity gives \(\det J=0\), while the bordered unit forces
generic pencil corank one and hence \(\operatorname{rank}J=1\). A primitive
generator \((P,Q)\) of \(\ker J\) gives the pencil-independent Hessian
kernel

\[
v=(0,0,P,Q)
\]

inside the fixed passive support plane. After exchanging active and passive
coordinate labels, this is precisely the two-component
quasitranslation-kernel stratum closed by HC4RSD5.

The genuinely new rank-two target consequently has
\(\operatorname{rank}\operatorname{Hess}_{z,w}B=1\). Its binary kernel
direction may still move rationally with both active and passive variables;
no univariate normal form is asserted here. This rank-one binary block,
together with the bordered-unit, Hessian integrability, and Piola equations,
is the next finite classification target. It is subsequently integrated and
closed by HC4RSD9 in
[HC4_QUADRATIC_RANK_TWO_PIVOTS.md](HC4_QUADRATIC_RANK_TWO_PIVOTS.md).

The rank-one pivot stratum is subsequently classified by HC4RSD10 in
[HC4_QUADRATIC_RANK_ONE_PIVOTS.md](HC4_QUADRATIC_RANK_ONE_PIVOTS.md).
Together these results close the full quadratic scalar-pivot branch with an
identically singular reduced Hessian pencil.

## 7. Reproduction

Run:

~~~bash
.venv/bin/python scripts/verify_hc4_quadratic_pivot_rank_obstruction.py
~~~

The command checks (2.3), (3.4), the odd-degree input, and the sharp
rank-two calibration. The degree-parity conclusion in (3.5) is the written
polynomial-degree argument above.
