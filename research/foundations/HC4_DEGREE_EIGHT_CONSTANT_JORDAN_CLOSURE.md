# Degree-eight scalar closure and constant-Jordan classification for `HC4`

## Status and scope

This note continues the synchronized scalar reverse-Schur packet of
[`HC4RSD20`](HC4_HIGHER_DEGREE_PENCIL_OBSTRUCTIONS.md) and the nilpotent
relative-pencil formulation of `HC4RSD17`.

It proves two exact results.

> **Theorem HC4RSD41 — complete scalar degree-eight closure.**  In the
> synchronized scalar reverse-Schur packet, every border coefficient whose
> leading binary form has degree at most eight has a fixed ruling.  Hence the
> packet reduces to `HC2` or to the exact `JC2` cotangent endpoint.

> **Theorem HC4RSD42 — constant nilpotent frame classification.**  Suppose
> the relative endomorphism
> \[
> N=S^{-1}T,
> \qquad S=\operatorname{Hess}\psi,
> \qquad T=\operatorname{Hess}A,
> \]
> is a constant nilpotent matrix and satisfies
> \(N^{\mathsf T}S=SN\).  Every nonzero Jordan type is either a triangular
> polynomial-automorphism pencil, an `HC2` endpoint, or exactly the `JC2`
> cotangent lift.  No constant nilpotent frame supplies an `HC4`-specific
> counterexample.

These statements do **not** prove unrestricted `HC4`.  Polynomially moving
nilpotent flags, direct four-variable packets outside reverse Schur, and the
`JC2` cotangent endpoint remain separate problems.

## 1. The scalar equation

After tangent-ruling synchronization, write the first nonbinary weighted
face as

\[
 c=f(x,y)+zg(x,y)+\frac{z^2}{2}q(x,y)+\cdots,
 \qquad \deg f=d.
\]

Its constant-in-\(z\) equation is

\[
 q\det B_f=b_g^{\mathsf T}\operatorname{adj}(B_f)b_g,
 \tag{1.1}
\]

where

\[
 B_f=
 \begin{pmatrix}
 0&f_x&f_y\\
 f_x&f_{xx}&f_{xy}\\
 f_y&f_{xy}&f_{yy}
 \end{pmatrix},
 \qquad
 b_g=(g,g_x,g_y)^{\mathsf T}.
\]

At a root of multiplicity \(m\), the valuation calculation in `HC4RSD27`
gives the resonance scalar

\[
 C_{d,m,e,n}
 =d^2m+d^2n^2-2demn-2dem-dm^2-dm
  +e^2m^2+2em^2+m^2.                                \tag{1.2}
\]

For \(d=8\), every proper root and every \(e\le6\) has no resonance below
\(\lceil m/2\rceil\).  Thus each root divides \(g\) to at least that order.
This converts the degree-eight problem into a finite root-partition census.

## 2. Complete same-weight faces

For \(e=6\), the complete weight-eight potential is

\[
 f+zg+\frac{z^2}{2}q+\frac{z^3}{6}r_2+\frac{z^4}{24}s. \tag{2.1}
\]

Keeping the \(z^3\) and \(z^4\) terms is essential.  Omitting either creates
false exceptional components in the `(4,2,2)` and four-double-root charts.

The exact partition calculation has the following structure.

1. Squarefree and one-double-root tops are already closed by
   `HC4RSD25--26`.
2. Every two- and three-root chart either has zero transverse radical in
   (1.1), or reaches a complete weighted face whose radical is the fixed
   cylinder.
3. Every weight-six partition has
   \[
   f=A^2B,\qquad g=AB,
   \qquad \deg A=2,\quad\deg B=4.
   \]
   Equation (1.1) forces
   \[
   q=\frac{31}{56}B
   \]
   and then
   \[
   A^2B''-3AA'B'+2AA''B+2(A')^2B=0.                \tag{2.2}
   \]
   The two affine quadratic normal forms for \(A\) leave no squarefree
   quartic \(B\).
4. For `(2,2,2,1,1)`, the values forced at the roots of the coprime cubic and
   quadratic factors are incompatible.
5. For four double roots, the only constant-Schur cross-ratios are
   \[
   \lambda\in\{-1,2,1/2\},
   \]
   one harmonic `PGL_2` orbit.  On the representative \(\lambda=-1\), the
   complete ideal is
   \[
   (a^3,a^2b,ab^2,b^3,r_0,r_1,r_2,s).
   \]
6. In the pure eighth-power chart, successive square coefficients remove
   every dependence on the third linear direction.

This exhausts all 22 partitions of eight and proves `HC4RSD41`.

## 3. Constant nilpotent relative frames

The constant nilpotent Jordan partitions of four are

\[
 [4],\qquad[3,1],\qquad[2,2],\qquad[2,1,1].
\]

### 3.1 Type `[4]`

The symmetry equation forces

\[
 S=
 \begin{pmatrix}
 0&0&0&a\\
 0&0&a&b\\
 0&a&b&c\\
 a&b&c&d
 \end{pmatrix},
 \qquad \det S=a^4.
\]

Hessian integrability and the unit determinant give

\[
\begin{aligned}
 \psi={}&aXW+aYZ+YU(W)+\frac{Z^2}{2}U'(W)
          +ZR'(W)+V(W),\\
 A={}&aYW+\frac a2Z^2+ZU(W)+R(W),
\end{aligned}                                       \tag{3.1}
\]

with \(a\in K^*\).  Every pencil gradient is inverted successively in the
order \(W,Z,Y,X\).

### 3.2 Type `[3,1]`

The integrated normal form is

\[
\begin{aligned}
 \psi={}&aXZ+\frac a2Y^2+YP'(Z)
          +\frac e2W^2+WQ(Z)+R(Z),\\
 A={}&aYZ+P(Z),
\end{aligned}                                       \tag{3.2}
\]

where \(a,e\in K^*\).  Its determinant is \(-a^3e\), and the inverse order
is \(Z,W,Y,X\).

### 3.3 Types `[2,1,1]` and `[2,2]`

Type `[2,1,1]` gives

\[
 \psi=aXY+C(Y,Z,W),\qquad A=\frac a2Y^2,
\]

and

\[
 \det\operatorname{Hess}\psi
 =-a^2\det\operatorname{Hess}_{Z,W}C,
\]

which is the `HC2` endpoint over the parameter \(Y\).

Type `[2,2]` gives

\[
 \psi=XH_Y(Y,W)+ZH_W(Y,W)+R(Y,W),\qquad A=H(Y,W),
\]

with

\[
 \det\operatorname{Hess}\psi
 =\det(\operatorname{Hess}H)^2.                     \tag{3.3}
\]

This is exactly the cotangent lift of the plane Keller map \(\nabla H\).
Together these four rows prove `HC4RSD42`.

## 4. Reproduction

Run

```bash
.venv/bin/python scripts/verify_hc4_degree_eight_partition_census.py
.venv/bin/python scripts/verify_hc4_degree_eight_four_five_root.py
.venv/bin/python scripts/verify_hc4_degree_eight_four_double.py
.venv/bin/python scripts/verify_hc4_degree_eight_pure_power.py
.venv/bin/python scripts/verify_hc4_constant_nilpotent_jordan.py
```

The checkers write their corresponding JSON records under
`artifacts/generated-results/` and verify every determinant, root-partition
ideal, harmonic orbit, pure-power square coefficient, and constant-Jordan
normal form used above.
