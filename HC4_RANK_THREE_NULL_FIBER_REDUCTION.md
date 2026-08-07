# Null affine-fiber reduction for the final rank-three `[4]` HC4 stratum

## Status and scope

Continue `HC4RSD64--65`.  Let

\[
S=\operatorname{Hess}\psi,
\qquad
T=\operatorname{Hess}A,
\qquad
\det(S+sT)=\delta\in K^*,
\]

and assume the relative nilpotent `N=S^{-1}T` has Jordan type `[4]`.
Thus `rank T=3` generically and `det T=0`.

> **Theorem HC4RSD66 — null affine-line fiber reduction.**
> Let `R` be a relation of minimum degree among the polynomial relations of
> `grad A`, and put
>
> \[
> k=(\nabla R)(\nabla A).
> \tag{0.1}
> \]
>
> Then `k` is nonzero and
>
> \[
> Tk=0,
> \qquad
> (k\cdot\nabla)k=0.
> \tag{0.2}
> \]
>
> The constant-Hessian pencil supplies the additional null identity
>
> \[
> k^{\mathsf T}Sk=0.
> \tag{0.3}
> \]
>
> Consequently, for `D=k dot grad`,
>
> \[
> D(\nabla A)=0,
> \qquad
> D^2\psi=0,
> \qquad
> D(D\psi)=0.
> \tag{0.4}
>
> On the generic rank-three locus, the fibers of `grad A` are exactly the
> affine lines
>
> \[
> x+t k(x).
> \tag{0.5}
>
> Thus every surviving `[4]` packet is a one-dimensional affine-line quotient
> on which `grad A` is constant and `psi` is affine.

> **Theorem HC4RSD67 — commuting inverse-Hessian frame.**
> Put
>
> \[
> G_s=(S+sT)^{-1}
> =G_0+sG_1+s^2G_2+s^3G_3.
> \tag{0.6}
> \]
>
> For every `s` and every standard basis vectors `e_i,e_j`, the polynomial
> vector fields
>
> \[
> X_i(s)=G_s e_i
> \]
>
> commute:
>
> \[
> [X_i(s),X_j(s)]=0.
> \tag{0.7}
> \]
>
> Therefore the cofactor flag of `HC4RSD64` satisfies, coefficient by
> coefficient,
>
> \[
> \sum_{a+b=m}[G_a e_i,G_b e_j]=0,
> \qquad 0\le m\le6.
> \tag{0.8}
> \]
>
> In particular the top rank-one line and all its interactions with the lower
> flag are constrained by a complete Lie-bracket hierarchy, strictly stronger
> than the row-divergence Piola identities alone.

The remaining `[4]` obstruction is therefore not an arbitrary moving singular
Hessian.  It is a null affine-line fibration equipped with a degree-three
commuting inverse-Hessian frame.

## 1. Associated quasi-translation

Since `rank Hess A=3`, the four components of `grad A` are algebraically
dependent only to transcendence degree three.  Hence a nonzero polynomial
relation

\[
R(\nabla A)=0
\tag{1.1}
\]

exists.  Differentiating gives

\[
((\nabla R)(\nabla A))^{\mathsf T}T=0.
\]

By symmetry of `T`, the vector `k` in (0.1) satisfies

\[
Tk=0.
\tag{1.2}
\]

For a relation of minimum degree, the associated vector is nonzero.  The
standard singular-Hessian/quasi-translation lemma of de Bondt gives

\[
(\mathcal Jk)k=0,
\tag{1.3}
\]

which is the second equation of (0.2).  Equivalently `x -> x+t k(x)` is a
polynomial quasi-translation and

\[
\nabla A(x+t k(x))=\nabla A(x).
\tag{1.4}
\]

## 2. The constant determinant makes the orbit null

For arbitrary square matrices `S,T`, the coefficient of `s^{n-1}` in
`det(S+sT)` is

\[
\operatorname{tr}(\operatorname{adj}T\,S).
\]

Here `n=4`, `rank T=3`, and `adj T` has rank one with image `ker T`.  On a
dense open set we may write

\[
\operatorname{adj}T=\rho k k^{\mathsf T}
\tag{2.1}
\]

for a nonzero rational scalar `rho`; replacing `k` by the associated
polynomial kernel vector does not change the line.

Since `det(S+sT)` is independent of `s`,

\[
0=[s^3]\det(S+sT)
 =\rho k^{\mathsf T}Sk.
\]

Thus

\[
k^{\mathsf T}Sk=0.
\tag{2.2}
\]

Because `Dk=0`,

\[
D^2\psi
=D(k\cdot\nabla\psi)
=(Dk)\cdot\nabla\psi+k^{\mathsf T}Sk
=0.
\tag{2.3}
\]

Likewise

\[
D(\nabla A)=Tk=0.
\tag{2.4}
\]

This proves (0.3)--(0.4).

## 3. Generic fibers are the orbit lines

The rank theorem gives generic fiber dimension

\[
4-\operatorname{rank}T=1
\]

for `grad A`.  Equation (1.4) places the entire affine line

\[
L_x=\{x+t k(x):t\in K\}
\]

inside the fiber through a generic point `x`.  A one-dimensional generic fiber
cannot contain a different positive-dimensional component through `x`.
Hence the generic fiber component is exactly `L_x`; after passing to the
function field of the image, this is the generic fiber itself.

The extra HC4 datum `psi` restricts to

\[
\psi(x+t k(x))=\psi(x)+tD\psi(x)
\tag{3.1}
\]

because `D^2 psi=0`.  Its slope `D psi` is itself constant along the orbit.

## 4. The commuting inverse-Hessian frame

For fixed `s`, define

\[
F_s=\psi+sA,
\qquad y=\nabla F_s(x).
\]

Since

\[
\det\operatorname{Hess}F_s=\delta\ne0,
\]

the gradient map is etale formally/locally, and

\[
\frac{\partial x}{\partial y}=G_s=(S+sT)^{-1}.
\]

Thus the columns

\[
X_i(s)=\sum_a(G_s)_{ai}\partial_{x_a}
\]

are precisely the pullbacks of the coordinate vector fields
`partial/partial y_i`.  Coordinate fields commute, proving (0.7).

Nilpotence of `N` gives

\[
G_s=(I-sN+s^2N^2-s^3N^3)S^{-1},
\tag{4.1}
\]

so the dependence on `s` has degree three.  Expanding the bracket in (0.7)
gives (0.8).  The top equations are

\[
[G_3e_i,G_3e_j]=0,
\tag{4.2}
\]

\[
[G_3e_i,G_2e_j]+[G_2e_i,G_3e_j]=0,
\tag{4.3}
\]

and

\[
[G_3e_i,G_1e_j]+[G_2e_i,G_2e_j]+[G_1e_i,G_3e_j]=0.
\tag{4.4}
\]

Since `im G_j=im N^j`, these are intrinsic differential constraints on the
nested Jordan/cofactor flag.

## 5. Proof boundary and external input

The only external ingredient in HC4RSD66 is the standard association between
a minimum-degree polynomial relation of a singular Hessian gradient and a
nonzero quasi-translation.  See M. de Bondt, *Quasi-translations and singular
Hessians*, arXiv:1501.05168, and the later small-rank Hessian work
arXiv:1609.03904.

Everything else is an exact determinant, Hessian, or coordinate-vector-field
identity.

The important negative boundary is also worth recording: pointwise `[4]`
nilpotence plus Hessian third-derivative symmetry does **not** force the
relative operator `N` to be Nijenhuis or Haantjes.  Hence nilpotent-Nijenhuis
normal forms cannot be imported without proving additional global identities.
