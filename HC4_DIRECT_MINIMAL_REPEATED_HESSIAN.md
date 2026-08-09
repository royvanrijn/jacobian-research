# Minimal repeated-Hessian obstruction for direct `HC4`

## Scope

Continue `HC4-DIR1--2`.  Let the top homogeneous form of a hypothetical
four-variable constant-Hessian polynomial have generic Hessian rank three.
After Gordan--Noether, write

\[
\Psi_D=f(x,y,z),\qquad
A=\operatorname{Hess}f,\qquad
\Delta=\det A\ne0.
\]

Assume every irreducible component of `Delta=0` has generic Hessian corank one.
Write

\[
\Delta=u\prod_\pi \pi^{m_\pi}
\]

and define the square-root ramification divisor

\[
R=\prod_\pi \pi^{\lfloor m_\pi/2\rfloor}.
\tag{0.1}
\]

> **Theorem HC4-DIR3 — ramification budget and minimal repeated case.**
> Let `j` be the first homogeneous descent order at which the fixed top kernel
> rotates.  Then
>
> \[
> j\le \deg R.
> \tag{0.2}
> \]
>
> In particular, if `deg R=1`, then `j=1`.  In that case the top ternary form
> `f` is either singular as a projective plane curve or is Thom--Sebastiani
> after a constant linear change.  Consequently a smooth non-TS rank-three
> top form with `deg R<=1` cannot occur in `HC4`.

The squarefree theorem `HC4-DIR2` is the case `deg R=0`.

---

## 1. Local valuation lemma

Let `pi` be an irreducible factor of `Delta` with multiplicity `m`, and work
in the DVR at the generic point of `pi=0`.  Generic corank one means that,
after invertible row/column operations over the DVR, the symmetric matrix `A`
has elementary valuations

\[
(0,0,m).
\]

Suppose

\[
\Delta\mid b^T\operatorname{adj}(A)b.
\tag{1.1}
\]

In a diagonal DVR model `diag(unit,unit,unit*pi^m)`, the only term of valuation
strictly below `m` is the square of the kernel component of `b`.  Therefore
that component has valuation at least

\[
\lceil m/2\rceil.
\]

It follows that

\[
A^{-1}b
\]

has a pole along `pi` of order at most

\[
\lfloor m/2\rfloor.
\tag{1.2}
\]

Squarefree components (`m=1`) contribute no pole.

Applying this at every component gives

\[
R A^{-1}b\in K[x,y,z]^3.
\tag{1.3}
\]

---

## 2. Homogeneity gives the ramification budget

Use the scaled block matrix

\[
M(t)=
\begin{pmatrix}
A(t)&b(t)\\b(t)^T&c(t)
\end{pmatrix},
\qquad
\det M(t)=\delta t^{4(D-2)}.
\]

Let

\[
b(t)=t^j b_j+O(t^{j+1}),\qquad b_j\ne0.
\]

As in `HC4-DIR2`, all lower coefficients of `c(t)` vanish through order
`2j-1`, and at order `2j`

\[
\Delta c_{2j}=b_j^T\operatorname{adj}(A)b_j.
\tag{2.1}
\]

The layer is affine in the top-kernel variable, so

\[
b_j=\nabla a_j
\]

for a homogeneous ternary form `a_j`.  Put

\[
X=A^{-1}b_j.
\]

Then `X` is homogeneous of degree `-j`.  Equation (1.3) says `RX` is a
polynomial homogeneous vector of degree

\[
\deg R-j.
\]

Since `X` is nonzero,

\[
\boxed{j\le\deg R}.
\]

This proves (0.2).

---

## 3. The minimal repeated case `deg R=1`

Now assume `deg R=1`.  Then, up to a scalar,

\[
R=\ell
\]

for a linear form `ell`, and the budget forces `j=1`.

Since `ell X` is a polynomial homogeneous vector of degree zero, there is a
nonzero constant vector `c` such that

\[
X=\frac{c}{\ell}.
\tag{3.1}
\]

Thus

\[
\nabla a=A X
=\frac{1}{\ell}\,A c
=\frac{1}{\ell}\nabla(D_cf).
\tag{3.2}
\]

Equivalently,

\[
\nabla(D_cf)=\ell\nabla a.
\tag{3.3}
\]

Take coordinates with `ell=z`.  The right side of (3.3) must itself be a
gradient.  Comparing mixed derivatives involving `z` gives

\[
a_x=a_y=0.
\]

Because `a` is homogeneous of degree `D-2`,

\[
a=\alpha z^{D-2}.
\tag{3.4}
\]

Integrating (3.3),

\[
D_cf=\beta z^{D-1}
\tag{3.5}
\]

for a constant `beta`.

There are two cases.

### 3.1 `D_c ell != 0`

Choose a linear coordinate `u` along `c`.  Since `ell` changes along `u`,
integrating (3.5) gives

\[
f=\gamma\ell^D+g(v_1,v_2)
\tag{3.6}
\]

after a constant linear change, where `v_1,v_2` are invariants of `D_c`.
This is Thom--Sebastiani.

### 3.2 `D_c ell = 0`

Now `ell` is an invariant of `D_c`, and integration gives

\[
f=\beta u\ell^{D-1}+g(\ell,v).
\tag{3.7}
\]

The projective point

\[
[u:\ell:v]=[1:0:0]
\]

is singular for `f=0` whenever `D>2`: all three first derivatives vanish
there.  Thus this alternative cannot occur when the top plane curve is
smooth.

This proves HC4-DIR3.

---

## 4. Interpretation

The non-squarefree Hessian locus is not itself enough to imply
Thom--Sebastiani: smooth non-TS ternary forms with a double Hessian component
exist.  HC4 supplies the additional descent vector `b_j`, and that extra datum
is what creates rigidity.

The correct complexity parameter for the direct rank-three attack is therefore
not degree and not merely reducibility of the Hessian curve.  It is

\[
\boxed{\sigma(f):=\deg R
=\deg\prod_\pi \pi^{\lfloor v_\pi(\Delta_f)/2\rfloor}.}
\]

A hypothetical HC4 counterexample must spend this ramification budget to
rotate its top kernel.  The first rotation order satisfies `j<=sigma(f)`.

The next target is `sigma(f)>=2`: classify the low-degree polynomial vector

\[
R(\operatorname{Hess}f)^{-1}\nabla a_j
\]

of degree `sigma-j` and use its curl/Hessian-integrability equations.  This
keeps the attack all-degree while increasing only the **ramification
complexity**, not the total degree.
