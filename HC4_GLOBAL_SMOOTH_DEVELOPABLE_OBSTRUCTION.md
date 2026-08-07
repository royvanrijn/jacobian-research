# Global smooth-developable obstruction for the scalar `HC4` packet

## Status and scope

This note replaces the positive-transverse-excess interpolation programme by a
global geometric obstruction.  It uses one irreducible smooth fiber and the
classical projective classification of developable surfaces.

> **Theorem HC4RSD52 — smooth polynomial developables are cylinders.**
> Let \(K\) be a characteristic-zero field and let
> \(c\in K[x,y,z]\).  Assume, after scalar extension to an algebraic closure,
>
> \[
> V(c_x,c_y,c_z)=\varnothing
> \tag{0.1}
> \]
>
> and
>
> \[
> \mathcal U(c):=
> (\nabla c)^{\mathsf T}
> \operatorname{adj}(\operatorname{Hess}c)\nabla c=0.
> \tag{0.2}
> \]
>
> Then there is a nonzero constant vector \(v\in K^3\) such that
>
> \[
> D_vc=v\cdot\nabla c=0.
> \tag{0.3}
> \]
>
> Equivalently, after a constant linear change of coordinates,
> \(c\in K[u_1,u_2]\): every level surface of \(c\) is a cylinder with one
> common affine ruling direction.

> **Corollary HC4RSD53 — complete scalar reverse-Schur closure.**
> Let
>
> \[
> \Psi(x,y,z,w)=w\,c(x,y,z)+D(x,y,z)
> \tag{0.4}
> \]
>
> have nonzero constant Hessian determinant.  Then \(c\) has a fixed
> direction.  Consequently every residual scalar reverse-Schur packet of
> `HC4RSD20`, in every degree and at every transverse excess, reduces to
> `HC2` or to the exact `JC2` cotangent endpoint.

This closes the **entire synchronized scalar reverse-Schur branch**, not merely
`h=0` or `h=1`.  It does not prove unrestricted `HC4`: non-scalar/coisotropic
pivots, polynomially moving matrix flags, direct four-variable constructions
outside this descent, and the `JC2` endpoint remain open.

## 1. Exact algebraic bridge from four to three variables

Put

\[
p=\nabla c,\qquad
E=\operatorname{Hess}D,\qquad
G=\operatorname{Hess}c.
\]

Then

\[
\operatorname{Hess}\Psi=
\begin{pmatrix}
E+wG&p\\
p^{\mathsf T}&0
\end{pmatrix}.
\tag{1.1}
\]

For every \(3\times3\) matrix \(M\),

\[
\det
\begin{pmatrix}
M&p\\
p^{\mathsf T}&0
\end{pmatrix}
=-p^{\mathsf T}\operatorname{adj}(M)p.
\tag{1.2}
\]

Therefore

\[
\det\operatorname{Hess}\Psi
=-p^{\mathsf T}\operatorname{adj}(E+wG)p.
\tag{1.3}
\]

The right side has degree at most two in \(w\), and

\[
[w^2]\det\operatorname{Hess}\Psi
=-p^{\mathsf T}\operatorname{adj}(G)p
=-\mathcal U(c).
\tag{1.4}
\]

If \(p\) vanished at any geometric point, the last row and column of (1.1)
would vanish there, contradicting the nonzero constant determinant.  Thus
(0.1) holds.  Constancy in \(w\) gives (0.2).

The rest of the argument is independent of the reverse-Schur construction.

## 2. A polynomial submersion is non-composite

Base-change to \(\overline K\).  Suppose

\[
c=P(g),\qquad \deg P>1.
\]

Choose a root \(\alpha\) of \(P'\).  Since \(g-\alpha\) is a nonconstant
polynomial, it has a zero over the algebraically closed field.  At any such
point,

\[
\nabla c=P'(g)\nabla g=0,
\]

contradicting (0.1).  Hence \(c\) is non-composite.

The classical Bertini--Krull irreducibility theorem now supplies a scalar
\(\lambda\in\overline K\) for which

\[
F=c-\lambda
\tag{2.1}
\]

is irreducible.  In fact the generic fiber of every non-composite polynomial
is irreducible.  Since \(\nabla F=\nabla c\) never vanishes, the affine surface

\[
S=V(F)\subset\mathbb A^3_{\overline K}
\tag{2.2}
\]

is smooth.

Only this one fiber will be used.

## 3. The fiber is developable

The bordered Hessian of \(F\) is

\[
\det
\begin{pmatrix}
\operatorname{Hess}F&\nabla F\\
(\nabla F)^{\mathsf T}&0
\end{pmatrix}
=-\mathcal U(c)=0.
\tag{3.1}
\]

The classical bordered-Hessian criterion says exactly that a smooth implicit
surface satisfying (3.1) is developable.  Equivalently, its projective tangent
plane map has rank at most one.

Let \(X\subset\mathbb P^3\) be the projective closure of \(S\).  Irreducibility
of \(F\) implies irreducibility of \(X\), and developability on the dense
affine open makes \(X\) a projective developable surface.

## 4. Smooth affine developables can only be cones at infinity

The classical projective classification gives three possibilities.

### 4.1 Plane

If \(X\) is a plane, the conclusion is immediate.

### 4.2 Non-conical tangent developable

Every irreducible developable surface in \(\mathbb P^3\) which is not a cone
is the tangential surface of its edge-of-regression curve \(E(X)\).  Unless
the surface is a plane, a tangent developable is singular along a dense open
part of that curve.

The affine surface \(S\) is smooth.  Hence the entire edge \(E(X)\) must lie
in the hyperplane at infinity \(H_\infty\).  But every projective tangent line
to a curve contained in the linear hyperplane \(H_\infty\) is itself
contained in \(H_\infty\).  Its tangential surface would therefore be
contained in \(H_\infty\), contradicting the fact that \(S\) is a nonempty
affine open of \(X\).

Thus the non-conical tangent case is impossible.

### 4.3 Cone

The remaining possibility is a cone.  A non-plane projective cone is singular
at its vertex.  Since \(S\) is smooth, the vertex must lie at infinity:

\[
p_\infty=[0:v_1:v_2:v_3],\qquad v\ne0.
\tag{4.1}
\]

For every \(x\in S\), the projective line joining \([1:x]\) to \(p_\infty\)
is contained in \(X\).  Its affine part is

\[
x+t v,\qquad t\in\overline K.
\]

Hence

\[
x+t v\in S
\quad\text{for all }x\in S,\ t\in\overline K.
\tag{4.2}
\]

So \(S\) is an affine cylinder with constant direction \(v\).

## 5. One cylindrical fiber forces the polynomial to be cylindrical

Differentiate (4.2) at \(t=0\).  The polynomial \(D_vc\) vanishes on the
irreducible hypersurface \(S=V(c-\lambda)\).  Therefore

\[
c-\lambda\mid D_vc.
\tag{5.1}
\]

But

\[
\deg D_vc<\deg(c-\lambda)
\]

unless \(D_vc=0\).  Thus

\[
D_vc=0
\]

over \(\overline K\).

Finally, the map

\[
K^3\longrightarrow K[x,y,z],
\qquad
v\longmapsto v\cdot\nabla c
\]

is \(K\)-linear.  Its kernel becomes nonzero after extension to
\(\overline K\), so its kernel was already nonzero over \(K\).  This proves
HC4RSD52.

## 6. Consequence for `HC4`

Section 1 applies HC4RSD52 directly to every polynomial of the residual form
(0.4).  Thus \(c\) has a fixed ruling.  The fixed-ruling theorem
`HC4RSD20` then gives precisely two endpoints:

1. an `HC2` packet; or
2. the cotangent lift of a plane Keller map, i.e. the exact `JC2` endpoint.

This proves HC4RSD53.

The degree-by-degree and transverse-excess calculations remain valuable as
independent elementary certificates, but they are no longer required for the
scalar branch.  In particular:

\[
\boxed{\text{there is no live scalar }h=1\text{ interpolation problem.}}
\]

Nor is there a live scalar \(h\ge2\) problem: the obstruction is global and
degree-free.

## 7. External inputs and proof boundary

The exact algebra in Section 1 is verified by the companion checker.  The
geometric proof uses the following classical inputs.

1. **Bertini irreducibility:** the generic fiber of a non-composite polynomial
   is irreducible.
2. **Bordered-Hessian criterion:** an implicit surface is developable exactly
   when its bordered Hessian vanishes on the surface.
3. **Projective developable classification:** every irreducible projective
   developable surface in \(\mathbb P^3\) is ruled; a non-cone is the
   tangential surface of its edge of regression.
4. **Tangent-edge singularity:** a non-planar tangent developable is singular
   along its generating curve.

Useful modern references are:

- J. Volčič, *Free Bertini's theorem and applications*, Proc. Amer. Math.
  Soc. 148 (2020), 3661--3671, DOI `10.1090/proc/15071`.
- S. Pérez-Díaz and L.-Y. Shen, *Determination and (re)parametrization of
  rational developable surfaces*, arXiv:`1305.2463`, especially Theorem 3.1.
- K. Kohn, B. Sturmfels and M. Trager, *Changing Views on Curves and
  Surfaces*, Acta Math. Vietnam. 43 (2018), 1--29, especially Theorem 2.2
  and the edge-of-regression discussion.
- G. Ishikawa, *Singularities of Developable Surfaces*, Hokkaido University
  Preprint Series 448 (1999).

## 8. Verification

Run

```bash
.venv/bin/python \
  scripts/verify_hc4_global_smooth_developable_obstruction.py
```

The checker verifies

\[
\det
\begin{pmatrix}
E+wG&p\\p^{\mathsf T}&0
\end{pmatrix}
=-p^{\mathsf T}\operatorname{adj}(E+wG)p,
\]

its \(w^2\)-coefficient, the gradient-zero contradiction, and the equality of
the bordered-Hessian determinant with \(-\mathcal U(c)\).  The Bertini and
projective-classification steps are external geometric theorems and are
declared as such rather than being represented by a CAS calculation.
