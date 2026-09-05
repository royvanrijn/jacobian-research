# Two-component quasi-translation kernels for HC4

## Status

This note continues HC4RSD4 from a fixed shear direction to a fixed
primitive direction supported in a constant two-plane. The coefficients
are initially allowed to depend on all four reduced variables.

> **Theorem HC4RSD5 (two-component quasi-translation obstruction).**
> Let \(K\) have characteristic zero and suppose
>
> \[
> \Phi(x,t)=\frac{\lambda}{2}t^2+tA(x)+B(x),
> \qquad \det\operatorname{Hess}\Phi=c\in K^\times,
> \tag{0.1}
> \]
>
> while the reduced pencil
> \(M(s)=\operatorname{Hess}_x(B+sA)\) is identically singular. Suppose its
> generic kernel line is independent of \(s\) and has a unimodular generator
>
> \[
> v=(P(x,y,z,w),Q(x,y,z,w),0,0)^{\mathsf T},
> \qquad (P,Q)=K[x,y,z,w].
> \tag{0.2}
> \]
>
> The Hessian Piola identity forces \(P,Q\in K[z,w]\). Then either the
> projective direction of \((P,Q)\) is constant, in which case HC4RSD1
> applies, or a constant linear change of the active coordinates transforms
> \(v\) into
>
> \[
> (L(z,w),1,0,0)^{\mathsf T}.
> \tag{0.3}
> \]
>
> The HC4RSD4 bordered-unit calculation then forces \(L\) to be a
> polynomial in one transverse linear form, and every scalar Schur
> descendant has an explicit triangular polynomial inverse. Hence every
> fixed primitive two-component kernel in a constant support plane has no
> collision.

The exact verifier is
[scripts/verify_hc4_two_component_quasitranslation_kernels.py](scripts/verify_hc4_two_component_quasitranslation_kernels.py),
and its generated ledger is
[artifacts/generated-results/hc4_two_component_quasitranslation_kernels.json](artifacts/generated-results/hc4_two_component_quasitranslation_kernels.json).

## 1. Piola eliminates active-variable dependence

The scalar Schur bordered unit makes the reduced pencil generically
corank one. With the primitive normalization (0.2), its adjugate is a
nonzero constant multiple of \(vv^{\mathsf T}\). The Hessian Piola identity

\[
\operatorname{div}(vv^{\mathsf T})=0
\tag{1.1}
\]

has two nonzero components:

\[
\begin{aligned}
2P P_x+Q P_y+P Q_y&=0,\\
P Q_x+Q P_x+2Q Q_y&=0.
\end{aligned}
\tag{1.2}
\]

If one of the two components is zero, unimodularity makes the other a unit
and the kernel direction is constant. Assume both are nonzero. Reduce the
first equation modulo \(P\). Since \((P,Q)=1\), it gives
\(P\mid P_y\), hence \(P_y=0\) by degree in \(y\). Similarly, the second
equation modulo \(Q\) gives \(Q_x=0\). The remaining equations are

\[
2P_x+Q_y=0,
\qquad
P_x+2Q_y=0.
\tag{1.3}
\]

Their coefficient determinant is \(3\), so characteristic zero gives

\[
P_x=P_y=Q_x=Q_y=0.
\tag{1.4}
\]

Thus \(P,Q\in K[z,w]\); transverse dependence is a consequence, not an
extra hypothesis.

## 2. Complete kernel integration

Write the active variables as \(x,y\) and the transverse variables as
\(u=(z,w)\). For any pencil member \(f=B+sA\), the kernel equation is

\[
\operatorname{Hess}(f)v=0,
\qquad
\delta=P(u)\partial_x+Q(u)\partial_y.
\tag{2.1}
\]

Choose \(U,V\in K[u]\) with \(UP+VQ=1\), and set

\[
\ell=Ux+Vy,
\qquad
r=Qx-Py.
\tag{2.2}
\]

Then \(\delta\ell=1\), \(\delta r=0\), and

\[
x=P\ell+Vr,
\qquad
y=Q\ell-Ur.
\tag{2.3}
\]

Every gradient component of \(f\) lies in the invariant ring \(K[u,r]\).
Put \(X=f_x\in K[u,r]\). Equality of the \(x,u_i\) mixed partials has an
\(\ell\)-coefficient

\[
(P Q_{u_i}-Q P_{u_i})X_r.
\tag{2.4}
\]

If the projective direction of \((P,Q)\) moves, (2.4) for \(z\) or \(w\)
forces \(X_r=0\); the same argument and \(f_{xy}=f_{yx}\) give
\((f_y)_r=0\). Exact integration now yields

\[
\boxed{
f=x a(z,w)+y b(z,w)+G(z,w),
\qquad
P\,da+Q\,db=0.
}
\tag{2.5}
\]

Conversely, (2.5) directly gives \(\operatorname{Hess}(f)v=0\). If the
projective direction is constant, unimodularity makes it a constant
kernel line, already covered by HC4RSD1.

## 3. The bordered unit produces a polynomial frame

Apply (2.5) to \(A\), writing its active coefficients as \(a_A,b_A\). The
rank-three bordered determinant from (0.1) makes the kernel-gradient
pairing a nonzero constant:

\[
P a_A+Q b_A=\alpha\in K^\times.
\tag{3.1}
\]

Differentiating (3.1) and using the kernel integral gives

\[
a_A\,dP+b_A\,dQ=0.
\tag{3.2}
\]

The row \((a_A,b_A)\) is nonzero, so (3.2) implies
\(dP\wedge dQ=0\). Choose a closed polynomial common composite \(H(z,w)\):

\[
P=P_0(H),
\qquad
Q=Q_0(H).
\tag{3.3}
\]

Because the projective direction moves, equations (3.1)--(3.2) uniquely
express \(a_A,b_A\) in \(K(H)\). Closedness of \(H\) gives
\(K[z,w]\cap K(H)=K[H]\), so write

\[
a_A=a_0(H),
\qquad
b_A=b_0(H).
\tag{3.4}
\]

In \(K[H]^2\), define

\[
r_0=(P_0,Q_0),
\qquad
r_1=(-b_0,a_0).
\tag{3.5}
\]

Equations (3.1)--(3.2) say

\[
\det\binom{r_0}{r_1}=\alpha,
\qquad
r_0'=S r_1
\tag{3.6}
\]

for some \(S\in K[H]\). Differentiating the constant determinant shows

\[
r_1'=T r_0
\tag{3.7}
\]

for some \(T\in K[H]\). Polynomiality of \(S,T\) follows directly from
the fact that the frame in (3.6) lies in
\(\operatorname{GL}_2(K[H])\).

## 4. Degree rigidity and reduction to HC4RSD4

Projective motion makes \(S\ne0\). If \(T\ne0\) as well, polynomial degree
in \(H\) gives

\[
\deg r_0-1=\deg S+\deg r_1,
\qquad
\deg r_1-1=\deg T+\deg r_0.
\tag{4.1}
\]

Adding the equations would give
\(-2=\deg S+\deg T\), a contradiction. Therefore \(T=0\), so \(r_1\) is
constant. For constants \(a_0,b_0,p_0,q_0\) and a nonconstant
\(L\in K[H]\), integration of (3.6) gives

\[
(P_0,Q_0)=(p_0,q_0)+L(H)(-b_0,a_0),
\qquad
a_0p_0+b_0q_0=\alpha.
\tag{4.2}
\]

Thus

\[
(P,Q)=(1,L(H))
\begin{pmatrix}
p_0&q_0\\
-b_0&a_0
\end{pmatrix}.
\tag{4.3}
\]

The constant matrix in (4.3) is invertible. A constant active-coordinate
change, followed if necessary by swapping the two active variables, turns
the kernel into \((L(H),1,0,0)\). This is precisely the HC4RSD4 input.
That theorem forces \(H\), and hence \(L(H)\), through one transverse
linear form and supplies the triangular inverse for every descendant.

## 5. Reproduction and frontier

Run:

~~~bash
.venv/bin/python scripts/verify_hc4_two_component_quasitranslation_kernels.py
# cleanup only: verify committed inputs and exact boundary
.venv/bin/python scripts/verify_hc4_two_component_quasitranslation_kernels.py --audit-existing-only
~~~

The command checks the Piola active-dependence gate, the complete Hessian
residual, the invariant-slice projective coefficient, the
constant-determinant polynomial frame, and the constant coordinate
reduction to the HC4RSD4 shear form.

The cleanup-only mode hash-checks the committed ledger and imported equation
helper without importing SymPy, replaying the identities, or rewriting the
artifact.

The theorem closes all fixed primitive kernel generators supported in a
constant two-plane. It does not classify fixed kernel generators with three
or four nonlinear components, or kernel lines that move with the pencil
parameter. Nonsingular reduced pencils with Schur-term cancellation and
matrix pivots with moving kernel planes are also outside its scope. The
nonzero-corner auxiliary constant-Hessian-pencil branch is consolidated
in `HC4MR1`, whose corrected maximal-motion closure uses HC4MRA1 and
[HC4MRA2](HC4_NEGATIVE_MOTION_POLYNOMIAL_OBSTRUCTION.md); nonlinear
zero-corner exact remainders and moving matrix planes remain open.
