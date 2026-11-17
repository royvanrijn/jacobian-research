# Ordinary-Laplacian GVC: low-rank extraction and Dvorsky-polarization obstructions

## 1. Status and scope

Work over a characteristic-zero field, normally \(\mathbb C\).  The target is a
polynomial \(P\) for the full-rank quadratic operator

\[
\Delta=\sum_i\partial_i^2
\]

such that

\[
\Delta^m(P^m)=0\quad(m\ge1),
\qquad
\Delta^m(P^{m+d})\ne0
\]

for one fixed \(d\ge1\) and infinitely many \(m\).  Full-rank quadratic
constant-coefficient operators are linearly equivalent over \(\mathbb C\), so
split formulas below are ordinary-Laplacian statements in disguised
coordinates.

This note does **not** give a global counterexample or prove ordinary GVC in
4--10 variables.  It closes four especially direct routes:

1. every affine-quadratic Hessian-nilpotent extraction from the five-variable
   Meng--Yang Schur family;
2. the complete naturally multigraded homogeneous-cubic seven-variable
   square-block completion of the Dvorsky seed;
3. the complete naturally multigraded homogeneous-cubic eight-variable
   split-pair completion;
4. the complete coarse homogeneous-cubic eight-variable two-square
   completion.

The exact checker is
`scripts/verify_ordinary_delta_frontier.py`.

## 2. The five-variable constant-Hessian family does not normalize to HN

Put \(u=1+x_1x_2\) and

\[
\begin{aligned}
A={}&y_1u^3+3x_1y_2u^2-x_1^3y_3,\\
B={}&y_1x_2^2u(4+3x_1x_2)
+y_2\bigl(x_2+3x_1x_2^2(4+3x_1x_2)\bigr)\\
&+y_3(2x_1-3x_1^2x_2).
\end{aligned}
\]

Consider the whole two-parameter Schur family

\[
\Psi_{\alpha,\beta}=A^2+\alpha A+\beta B.
\]

For \(\beta\ne0\), Schur descent gives

\[
\det\operatorname{Hess}\Psi_{\alpha,\beta}=8\beta^4.
\]

The published representative is \((\alpha,\beta)=(13,2)\), with determinant
\(128\).

### Theorem 2.1 — no affine-quadratic HN normalization

There are no invertible affine coordinates \(x=Sz+b\) and no constant
quadratic polynomial \(q(z)\) for which

\[
P(z)=\Psi_{\alpha,\beta}(Sz+b)-q(z)
\]

is Hessian nilpotent.  In fact, for every constant symmetric matrix \(K\),

\[
\operatorname{tr}\!\left(K\operatorname{Hess}\Psi_{\alpha,\beta}\right)
\in k
\quad\Longrightarrow\quad K=0.
\]

#### Proof

Write \(K=(k_{ij})_{1\le i,j\le5}\), symmetrically, in the variable order
\((x_1,x_2,y_1,y_2,y_3)\).  In the polynomial

\[
T_K=\operatorname{tr}(K\operatorname{Hess}\Psi_{\alpha,\beta}),
\]

the following coefficients are, in order,

\[
\begin{array}{c|c}
\text{monomial}&\text{coefficient}\\
\hline
x_1^6x_2^6&2k_{33}\\
x_1^6x_2^5y_1&24k_{23}\\
x_1^6x_2^5&12k_{34}\\
x_1^6x_2^4y_1^2&30k_{22}\\
x_1^6x_2^4y_1&60k_{24}\\
x_1^6x_2^4&18k_{44}\\
x_1^6x_2^3&-4k_{35}\\
x_1^6x_2^2y_1&-12k_{25}\\
x_1^6x_2^2&-12k_{45}\\
x_1^6&2k_{55}\\
x_1^5x_2^6y_1&24k_{13}\\
x_1^5x_2^5y_1^2&72k_{12}\\
x_1^5x_2^5y_1&72k_{14}\\
x_1^5x_2^3y_1&-24(k_{15}-10k_{24})\\
x_1^4x_2^6y_1^2&30k_{11}.
\end{array}
\]

Constancy of \(T_K\) therefore forces all fifteen entries of \(K\) to zero.

Now suppose \(P(z)=\Psi(Sz+b)-q(z)\) were HN.  Nilpotency of
\(\operatorname{Hess}P\) implies trace zero, hence

\[
\operatorname{tr}\left(SS^T\operatorname{Hess}\Psi(Sz+b)\right)
=\operatorname{tr}(\operatorname{Hess}q),
\]

a constant.  Invertible affine substitution preserves nonconstancy, so the
preceding result gives \(SS^T=0\), impossible for invertible \(S\).  \(\square\)

Thus constant Hessian determinant is not enough to extract an ordinary
Laplacian witness by subtracting a quadratic; the full Hessian pencil, not one
determinant value, is the missing condition.

## 3. Seven variables: complete square-block cubic lift

Start with the Dvorsky cubic

\[
P_0=(t+c)(ad+bt)
\]

and the full-rank quadratic operator

\[
\Delta_7=\partial_a\partial_d-
\partial_b\partial_c+
\partial_t\partial_s+
\partial_u^2.
\]

Use the integral multigrading

\[
\begin{aligned}
&\deg a=(2,0,0),\quad \deg d=(0,2,0),\\
&\deg b=\deg s=(2,2,-2),\\
&\deg c=\deg t=(0,0,2),\quad \deg u=(1,1,0).
\end{aligned}
\]

Every term of \(\Delta_7\) has degree \((-2,-2,0)\), and \(P_0\) has degree
\((2,2,2)\).  The complete cubic of this weight restricting to \(P_0\) on
\(s=u=0\) is

\[
F_7=P_0+q_0tu^2+q_1st^2+q_2cu^2+q_3cst+q_4c^2s.
\]

(The only further same-weight cubic is \(bc^2\), which would alter the required
hyperplane restriction.)

### Theorem 3.1 — moment-three obstruction

No coefficient choice satisfies

\[
\Delta_7^m(F_7^m)=0\qquad(m=1,2,3).
\]

#### Proof

The first moment gives

\[
q_0=-q_1,\qquad q_3=-2q_2-1.
\]

After this substitution, the second moment is

\[
\Delta_7^2(F_7^2)=4\bigl(E_1c^2-2E_2ct+E_3t^2\bigr),
\]

where

\[
\begin{aligned}
E_1&=6q_2^2+2q_2-2q_4+1,\\
E_2&=6q_1q_2+q_1-q_2+2q_4-1,\\
E_3&=6q_1^2+4q_2+3.
\end{aligned}
\]

Eliminating \(q_1,q_2\) from \((E_1,E_2,E_3)\) gives

\[
g(q_4)=72q_4^3+66q_4^2+20q_4-25=0.
\]

Modulo the same ideal, the normalized \(c^3\)-coefficient of the third moment
is

\[
[ c^3 ]\frac{\Delta_7^3(F_7^3)}{36}
\equiv-\frac23(2q_4-1)(4q_4+5).
\]

Put \(h=(2q_4-1)(4q_4+5)\).  The explicit Bézout identity

\[
(32q_4+34)g-(288q_4^2+354q_4+275)h=525
\]

shows that \(g=h=0\) is impossible in characteristic zero.  \(\square\)

## 4. Eight variables: complete split-pair cubic lift

Use

\[
\Delta_8=\partial_a\partial_d-
\partial_b\partial_c+
\partial_t\partial_s+
\partial_u\partial_v.
\]

Give \(a,u\) one weight, \(d,v\) its opposite half, \(b,s\) the common
return weight, and \(c,t\) the phase weight:

\[
\deg a=\deg u=(1,0,0),\quad
\deg d=\deg v=(0,1,0),
\]

\[
\deg b=\deg s=(1,1,-1),\quad
\deg c=\deg t=(0,0,1).
\]

The complete cubic of weight \((1,1,1)\), restricting to \(P_0\) when
\(s=u=v=0\), is

\[
\begin{aligned}
F_8=P_0
&+p_0tuv+p_1st^2+p_2dtu+p_3cuv+p_4cst\\
&+p_5cdu+p_6c^2s+p_7atv+p_8acv.
\end{aligned}
\]

### Theorem 4.1 — a unit at moment four

No coefficient choice satisfies

\[
\Delta_8^m(F_8^m)=0\qquad(m=1,2,3,4).
\]

#### Proof

The first moment gives

\[
p_0=-2p_1,\qquad p_3=-p_4-1.
\]

Rename

\[
x=p_1,\,r=p_2,\,q=p_4,\,h=p_5,\,z=p_6,\,\ell=p_7,\,w=p_8.
\]

Then

\[
\Delta_8^2(F_8^2)=4(E_1c^2+E_2ct+E_3t^2),
\]

\[
\Delta_8^3(F_8^3)=36(E_4c^3+E_5c^2t+E_6ct^2+E_7t^3),
\]

where

\[
\begin{aligned}
E_1={}&q^2+q+hw-2z+1,\\
E_2={}&4xq+2x+rw-q+h\ell-4z+1,\\
E_3={}&4x^2+r\ell-2q+1,\\
E_4={}&4xqz+2xz+rzw-q^2-qhw-3qz-q+hz\ell-4z^2+z,\\
E_5={}&8x^2z+4xq^2-2xq-2xhw-4xz-2x+2rz\ell\\
&-4q^2-12qz-q-hw+8z-1,\\
E_6={}&12x^2q-2x^2-xrw-7xq-xh\ell-12xz-x+rq\ell\\
&-8q^2+2q-2hw+12z-2,\\
E_7={}&8x^3-12xq+r\ell-rw+3q-h\ell+4z-1.
\end{aligned}
\]

Let

\[
\begin{aligned}
T={}&[t^4]\frac{\Delta_8^4(F_8^4)}{576}\\
={}&36x^4+7x^2r\ell-58x^2q+7x^2+xr\ell-xrw+19xq-xh\ell\\
&+20xz-x+r^2\ell^2-3rq\ell+2r\ell+11q^2-7q+hw-10z+2.
\end{aligned}
\]

A direct polynomial identity is

\[
T-1=\sum_{i=1}^7M_iE_i,
\]

with

\[
\begin{aligned}
M_1&=8x^2+4x-12q-1,\\
M_2&=-2x^2-2xq+3q+8z+1,\\
M_3&=3x^2+2xq-8xz-2x+r\ell+2q^2-2q-4z,\\
M_4&=-8,\\
M_5&=2(2x+1),\\
M_6&=-2(x+q+1),\\
M_7&=3x+3q+1.
\end{aligned}
\]

Thus moments two and three force \(E_1=\cdots=E_7=0\), but then \(T=1\), so
\([t^4]\Delta_8^4(F_8^4)=576\ne0\).  \(\square\)

## 5. Eight variables: complete coarse two-square lift

There is another complete eight-variable cubic architecture.  Give both
auxiliary variables \(z=(u,v)\) the half-weight \((1,1,0)\), and use

\[
\Delta_8'=\partial_a\partial_d-
\partial_b\partial_c+
\partial_t\partial_s+
\partial_u^2+
\partial_v^2.
\]

For arbitrary symmetric \(2\times2\) matrices \(A,B\), the complete cubic
with the required hyperplane restriction and vanishing first moment is

\[
\begin{aligned}
F_{A,B,\gamma}=P_0
&+s\bigl(-\operatorname{tr}(A)t^2
-(1+2\operatorname{tr}(B))ct+\gamma c^2\bigr)\\
&+t\,z^TAz+c\,z^TBz.
\end{aligned}
\]

### Theorem 5.1 — coarse moment-four obstruction

No \((A,B,\gamma)\) makes the first four pure moments vanish.

#### Exact certificate

The coefficient ideal of moments two and three contains

\[
g(\gamma)=8\gamma^3+12\gamma^2+10\gamma+5.
\]

Modulo that ideal, the five degree-four coefficients are nonzero rational
multiples of

\[
h(\gamma)=4\gamma^2+10\gamma+15
\]

(the middle coefficient carries an additional factor \(3\)).  Finally

\[
g+(2-2\gamma)h=35.
\]

Hence moments two and three force \(g=0\), while moment four would force
\(h=0\), an impossibility.  The checker reconstructs the complete coefficient
ideal from \(F_{A,B,\gamma}\); no random specialization is used.

## 6. Consequences for the next search

These closures separate the remaining plausible routes sharply.

* The five-variable constant-Hessian counterexample cannot be converted into
  an HN polynomial by affine normalization and quadratic subtraction.
* Adding the first one or two natural quadratic polarization blocks to the
  Dvorsky cubic does not work inside the complete compatible cubic spaces.
* A successful 4--10-variable construction must therefore leave at least one
  of these hypotheses: use non-multigraded corrections, degree at least four,
  a nonlinear polarization/specialization, a different seed, or a larger
  auxiliary architecture.

The most concrete next finite target is the **complete non-equivariant cubic
seven-variable lift**.  After fixing the Dvorsky hyperplane restriction it has
many more coefficients than the five-parameter square-block model, but the
present moment-three Bézout certificate supplies a strong elimination pivot.
The next structurally cleaner target is a nine-variable two-block lift, where
one should work in trace-word invariants rather than raw coefficients.

## 7. References used for orientation

* W. Zhao, *Hessian Nilpotent Polynomials and the Jacobian Conjecture*,
  arXiv:math/0409534.
* G. Meng and L. Yang, *A five-variable counterexample to the Hessian
  conjecture, and the low-dimensional status of the Jacobian and Hessian
  conjectures*, arXiv:2607.22198v2.
* Repository notes `DVORSKY_ONE_PAIR_SCHUR_OBSTRUCTION.md`,
  `TWO_VARIABLE_GVC_REPRESENTATION_PROGRAM.md`, and
  `THREE_VARIABLE_HOMOGENEOUS_GVC_COUNTEREXAMPLE.md`.
