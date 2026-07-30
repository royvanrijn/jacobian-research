# Universal relative quadratic-gauge Keller map

This note packages the supplied-polynomial realization into one relative
geometric object.  It also separates three notions which should not be
conflated:

1. a relative Keller map whose coefficients vary in a parameter scheme;
2. its marked finite-etale fiber incidence space;
3. a choice-free construction on the stack of abstract finite-etale
   algebras.

The first two are explicit.  The third is a descent problem, and the
quadratic-gauge family does not generically descend in stable
left--right moduli when `N>=5`.

Throughout, work over `Q` and fix `N>=3`.  Hasse coefficients are used in
the formulas, so the same construction works over any base on which `6` is
invertible after replacing the ordinary third derivative by the third Hasse
derivative.

## 1. The literal coefficient atlas

Let

\[
 {\cal P}(T)=T^N+c_{N-1}T^{N-1}+\cdots+c_0
\]

be the universal monic polynomial, let `a` be another variable, and write

\[
 d_j={\cal P}^{[j]}(a).
\]

Define

\[
 U_N=
 \operatorname{Spec}
 \mathbb Q[c_0,\ldots,c_{N-1},a,
   \Delta({\cal P})^{-1},d_1^{-1},d_3^{-1}].           \tag{1.1}
\]

On `U_N`, put

\[
 G(S)={\cal P}(a+S)-{\cal P}(a)
     =d_1S+d_2S^2+\cdots+d_NS^N.                       \tag{1.2}
\]

Use the quadratic-gauge formulas with this seed.  In source coordinates
`(x,y,z)`, set

\[
 t=1+xy,\qquad
 q=t^2z+\frac{d_1}{d_3}y^2(1+3t),                     \tag{1.3}
\]

and

\[
\begin{aligned}
 \Pi={}&tq,\\
 B={}&y+3\frac{d_3}{d_1}xq
       +2\frac{d_2}{d_1}tq
       +\sum_{j=4}^N j\frac{d_j}{d_1}
          t^2x^{j-2}q^j,\\
 C={}&x(5-3t)-\frac{d_3}{d_1}x^3z
       -\sum_{j=4}^N(j-2)\frac{d_j}{d_1}(xq)^j.
                                                               \tag{1.4}
\end{aligned}
\]

The normalized map is

\[
 \widetilde F_G=(\Pi,-B/2,C).
\]

Leaving the parameters fixed defines

\[
 \mathfrak F_N:
 U_N\times\mathbb A^3\longrightarrow
 U_N\times\mathbb A^3,\qquad
 (u,v)\longmapsto(u,\widetilde F_{G_u}(v)).             \tag{1.5}
\]

The vertical Jacobian is one by the quadratic-gauge theorem.  In the
coordinates of (1.1), the full Jacobian matrix has the block form

\[
 \begin{pmatrix}
 I_{N+1}&0\\
 *&D_v\widetilde F_G
 \end{pmatrix},
\]

and therefore

\[
 \boxed{\det D\mathfrak F_N=1.}                        \tag{1.6}
\]

The distinguished target section is

\[
 y_{{\cal P},a}
 =\left(1,0,-\frac{2{\cal P}(a)}{d_1}\right).          \tag{1.7}
\]

At this section the inverse equation is

\[
 G(S)-\frac{d_1}{2}
 \left(-\frac{2{\cal P}(a)}{d_1}\right)
 ={\cal P}(a+S).                                       \tag{1.8}
\]

Since the discriminant is inverted, derivative localization disappears.
Consequently there is a canonical cartesian fiber identification

\[
\boxed{
 (U_N\times\mathbb A^3)
 \mathop{\times}_{\mathfrak F_N,\,
 U_N\times\mathbb A^3,\,y}
 U_N
 \simeq
 \operatorname{Spec}_{U_N}
 \frac{{\cal O}_{U_N}[S]}{({\cal P}(a+S))}
 \simeq
 \operatorname{Spec}_{U_N}
 \frac{{\cal O}_{U_N}[T]}{({\cal P}(T))}.}             \tag{1.9}
\]

The second isomorphism is `T=a+S`.  It is an isomorphism of finite-etale
`U_N`-schemes of rank `N`, not only a bijection on geometric points.

### 1.1 The universal ordered-collision cover

The collision-algebra and divided-difference interface used here is credited
to Chloe van der Vlugt's *Collision Ideals and Off-Diagonal Sheets*.  The
[external audit and credit ledger](COLLISION_IDEALS_EXTERNAL_AUDIT.md)
records the separately attributed manuscript, public Lean repository,
license, AI-assistance disclosure, and exact formalization boundary.  The
universal ordered-root application below is a deduction in this repository,
not a claim that it appears in that paper.

Put

\[
 R_N={\cal O}(U_N),\qquad
 B_N=R_N[T]/({\cal P}(T)).
\]

By (1.9), the ordered self-collision fiber of the distinguished Keller fiber
has coordinate algebra

\[
\boxed{
 C_N=B_N\otimes_{R_N}B_N
 \simeq
 \frac{R_N[T_1,T_2]}
      {({\cal P}(T_1),{\cal P}(T_2))}.}                \tag{1.10}
\]

Diagonal restriction is multiplication

\[
 \mu_N:C_N\longrightarrow B_N,\qquad
 f(T_1,T_2)\longmapsto f(T,T),                        \tag{1.11}
\]

and the fiberwise collision obstruction is
`\mathfrak o_N=\ker(\mu_N)`.  Introduce the divided difference

\[
 D_{\cal P}(T_1,T_2)
 =\frac{{\cal P}(T_2)-{\cal P}(T_1)}{T_2-T_1}.        \tag{1.12}
\]

The discriminant is a unit on `U_N`, so the class of
`${\cal P}'(T_1)$` is a unit in `B_N`.  If `\eta(T_1)` denotes its inverse,
then

\[
 e_\Delta=\eta(T_1)D_{\cal P}(T_1,T_2)\in C_N         \tag{1.13}
\]

is the canonical diagonal idempotent: it is one modulo `T_2-T_1` and zero
modulo `D_{\cal P}`.  Since

\[
 {\cal P}(T_2)-{\cal P}(T_1)
 =(T_2-T_1)D_{\cal P}(T_1,T_2),
\]

the Chinese remainder theorem gives the exact splitting

\[
\boxed{
 C_N\simeq B_N\times B_N^{\mathrm{off}},\qquad
 B_N^{\mathrm{off}}
 =\frac{R_N[T_1,T_2]}
        {({\cal P}(T_1),D_{\cal P}(T_1,T_2))}.}        \tag{1.14}
\]

Under (1.14), `e_\Delta=(1,0)`, multiplication is projection to the first
factor, and

\[
 \mathfrak o_N=(1-e_\Delta)C_N\simeq B_N^{\mathrm{off}}.
                                                               \tag{1.15}
\]

All three algebras are finite etale over `R_N`, with locally free ranks

\[
\boxed{
 \operatorname{rk}C_N=N^2,\qquad
 \operatorname{rk}B_N=N,\qquad
 \operatorname{rk}B_N^{\mathrm{off}}=N(N-1).}         \tag{1.16}
\]

Functorially, (1.10) represents ordered pairs of points of the selected
fiber; (1.14) separates equal pairs from ordered distinct pairs.  Thus every
compiled full fiber `Spec(A)` of the absolute map `{\cal U}_N` has the
canonical collision decomposition

\[
 \operatorname{Spec}(A\otimes_kA)
 =\Delta(\operatorname{Spec}A)
 \sqcup\operatorname{Spec}(A^{\mathrm{off}}),          \tag{1.17}
\]

with off-diagonal rank `N(N-1)`.  Special fibers may split further.  The
generic connectedness, normal-closure boundary, and higher distinct-root
tower are determined by monodromy in Section 5.

## 2. Exact normalization of the proposed base

The base (1.1) has one redundant translation parameter.  Define

\[
 H(S)=\frac{{\cal P}(a+S)}{{\cal P}'(a)}
 =h_0+S+h_2S^2+h_3S^3+\cdots+h_NS^N.                  \tag{2.1}
\]

Here

\[
 h_j=\frac{d_j}{d_1},\qquad
 h_3h_N\ne0.
\]

Let

\[
 V_N=\operatorname{Spec}
 \mathbb Q[h_0,h_2,\ldots,h_N,
 h_3^{-1},h_N^{-1},\Delta(H)^{-1}].                    \tag{2.2}
\]

Then there is an exact isomorphism

\[
 \boxed{U_N\simeq V_N\times\mathbb A^1_a.}             \tag{2.3}
\]

The forward map is `({\cal P},a)\mapsto(H,a)`.  Its inverse is

\[
 {\cal P}(T)=h_N^{-1}H(T-a).                           \tag{2.4}
\]

Indeed (2.4) is monic and
`${\cal P}'(a)=h_N^{-1}$`, so applying (2.1) recovers `H`.
Thus `V_N`, of dimension `N`, is the normalized marked-presentation
space, while the literal `U_N` has dimension `N+1`.

## 3. Sharp `N-3`-parameter relative map

The constant, quadratic, and cubic coefficients of `H` can be placed in
the three target coordinates.  This gives a smaller universal relative
map.

For `N=3`, put `K_3=Spec(Q)`.  For `N>=4`, put

\[
 K_N=\operatorname{Spec}
 \mathbb Q[u_4,\ldots,u_N,u_N^{-1}],                   \tag{3.1}
\]

so

\[
 \dim K_N=N-3.
\]

Over `K_N`, use the seed

\[
 G_u(S)=S+S^3+\sum_{j=4}^Nu_jS^j.                     \tag{3.2}
\]

Equations (1.3)--(1.4), now with `d_1=d_3=1`, `d_2=0`,
and `d_j=u_j`, define a single relative Jacobian-one map

\[
 \mathcal K_N:
 K_N\times\mathbb A^3\longrightarrow
 K_N\times\mathbb A^3.                                \tag{3.3}
\]

Write the three normalized target coordinates as `(\pi,b,c)`, where the
second coordinate is `-B/2`.  The inverse equation of (3.3) is

\[
 E_{u,\pi,b,c}(S)
 =S+bS^2+\pi S^3+\sum_{j=4}^Nu_j\pi^jS^j-\frac c2.
                                                               \tag{3.4}
\]

Let `I_N` be the open subscheme of

\[
 K_N\times\mathbb G_{m,\pi}\times\mathbb A^2_{b,c}
\]

where the discriminant of (3.4) is nonzero.  Then

\[
 \mathcal K_N^{-1}(\pi,b,c)
 \simeq\operatorname{Spec}{\cal O}_{I_N}[S]/(E_{u,\pi,b,c})
                                                               \tag{3.5}
\]

over `I_N`.  The change of variables

\[
 h_0=-c/2,\qquad h_2=b,\qquad h_3=\pi,\qquad
 h_j=u_j\pi^j\quad(j\ge4)                              \tag{3.6}
\]

is an isomorphism

\[
 \boxed{I_N\simeq V_N.}                                \tag{3.7}
\]

Its inverse is

\[
 \pi=h_3,\qquad b=h_2,\qquad c=-2h_0,\qquad
 u_j=\frac{h_j}{h_3^j}.                                \tag{3.8}
\]

Therefore (3.3) is a universal relative Keller map over only `N-3`
parameters, and its moving targets recover every normalized squarefree
degree-`N` presentation.

For the original pair `({\cal P},a)`, the specialization is

\[
\begin{aligned}
 \pi&=\frac{{\cal P}^{[3]}(a)}{{\cal P}'(a)},\\
 b&=\frac{{\cal P}^{[2]}(a)}{{\cal P}'(a)},\\
 c&=-\frac{2{\cal P}(a)}{{\cal P}'(a)},\\
 u_j&=
 \frac{{\cal P}^{[j]}(a)\,{\cal P}'(a)^{j-1}}
      {{\cal P}^{[3]}(a)^j}\qquad(j\ge4).
                                                               \tag{3.9}
\end{aligned}
\]

Substitution in (3.4) gives

\[
 E_{u,\pi,b,c}(S)
 =\frac{{\cal P}(a+S)}{{\cal P}'(a)}.                  \tag{3.10}
\]

This is the requested compiler in family form.

The count `N-3` is sharp for dominance over normalized polynomial
presentations in the following precise sense.  If a family over an integral
scheme `X`, together with its three target coordinates, induces a dominant
rational map

\[
 X\times\mathbb A^3\dashrightarrow V_N,
\]

then

\[
 \dim X+3\ge\dim V_N=N,
 \qquad\text{hence}\qquad
 \boxed{\dim X\ge N-3.}                                \tag{3.11}
\]

This is a statement about dominance over *presentations*.  It is not a
lower bound for an atlas of `BS_N`: the stack `BS_N` has a zero-dimensional
etale atlas.

### 3.1 One absolute atomic map in every rank

The parameter promotion used below in degree five is not special to degree
five.  Let

\[
 {\bf u}=(u_4,\ldots,u_N)
\]

for `N>=4`, and let `{\bf u}` be the empty tuple for `N=3`.  Retain these
parameters as unchanged source and target coordinates and define

\[
\boxed{
 {\cal U}_N({\bf u},x,y,z)
 =({\bf u},{\cal K}_{N,{\bf u}}(x,y,z)):
 \mathbb A^N_{\mathbb Q}\longrightarrow\mathbb A^N_{\mathbb Q}.
}
\]

Here `\mathcal K_{N,{\bf u}}` is given by (1.3)--(1.4) with
`d_1=d_3=1`, `d_2=0`, and `d_j=u_j`.  Thus this is one polynomial map on
the whole affine space; no parameter is inverted in its definition.  The
full Jacobian is block triangular with an `(N-3)`-by-`(N-3)` identity block
and the vertical determinant-one block, so

\[
 \det D{\cal U}_N=1.
\]

At a target `({\bf u},\pi,b,c)` its generic inverse equation is

\[
\boxed{
 E_{{\bf u},\pi,b,c}(S)
 =S+bS^2+\pi S^3+\sum_{j=4}^Nu_j\pi^jS^j-\frac c2.
}
\]

Its leading coefficient is `u_N\pi^N` for `N>=4` and `\pi` for `N=3`.
The reconstruction in (3.5) therefore proves

\[
 \operatorname{gdeg}({\cal U}_N)=N.
\]

Let `k` be any characteristic-zero field and let `A` be a rank-`N` finite
etale `k`-algebra.  Since `k` is infinite, `A` has a primitive element, so

\[
 A\simeq k[T]/(P(T))
\]

for a monic squarefree polynomial `P` of degree `N`.  Choose `a\in k`
outside the finite zero set of
`P'(a)P^{[3]}(a)`.  This is possible because both factors are nonzero
polynomials.  Write

\[
 \frac{P(a+S)}{P'(a)}
 =h_0+S+h_2S^2+h_3S^3+\cdots+h_NS^N
\]

and take

\[
 \pi=h_3,\qquad b=h_2,\qquad c=-2h_0,\qquad
 u_j=\frac{h_j}{h_3^j}\quad(4\le j\le N).
\]

Then the displayed inverse equation is exactly `P(a+S)/P'(a)`, and (3.5)
identifies the **complete** target fiber with

\[
 \operatorname{Spec}k[S]/(P(a+S))
 \simeq\operatorname{Spec}(A).
\]

On the generic target, the birational coefficient change (3.6) identifies
this inverse equation with the universal normalized root polynomial over
`V_N`.  Thus promoting the parameters does not change the generic inverse
extension of the relative map, and Section 5 gives geometric monodromy
`S_N`.  Its natural degree-`N` action is primitive.  The
[primitive-monodromy atomicity theorem](PRIMITIVE_MONODROMY_ATOMICITY.md)
therefore makes `{\cal U}_N` atomic after every characteristic-zero
extension of constants and after every identity stabilization.

> **Universal absolute atomic-map theorem.**  For every `N>=3`, the single
> explicit `Q`-defined map
> \[
> {\cal U}_N:\mathbb A^N\longrightarrow\mathbb A^N
> \]
> above has determinant one, geometric degree `N`, and geometric monodromy
> `S_N`.  It is absolutely and stably atomic.  After base change to any
> characteristic-zero field `k`, its complete target fibers realize every
> rank-`N` finite etale `k`-algebra.

The theorem is optimal with respect to the full-fiber rank spectrum in the
following limited sense.  Rank one is represented by an automorphism, while
the [complete rank classification](FINITE_ETALE_KELLER_FIBERS.md) excludes
full rank two in characteristic zero.  It does not claim that ambient
dimension `N` is minimal.  It also makes no canonical choice of the target
attached to an abstract algebra: the primitive element and the translation
`a` remain presentation choices, exactly as in the atlas discussion of
Section 4.

The
[adversarial audit](UNIVERSAL_ATOMIC_MAP_ADVERSARIAL_AUDIT.md)
checks every implication separately, records the parameter and separability
failure boundaries, and supplies connected, split, and product witness cards.

### 3.2 Explicit quintic coordinates

For `N=5`, the relative construction extends across `u_5=0` to one
polynomial Keller self-map of affine five-space.  Write the source
coordinates as `(u,v,x,y,z)` and put

\[
\begin{aligned}
 t&=1+xy,\\
 q&=t^2z+y^2(1+3t),\\
 B&=y+3xq+4ut^2x^2q^4+5vt^2x^3q^5,\\
 C&=x(5-3t)-x^3z-2u(xq)^4-3v(xq)^5.
\end{aligned}
\]

Then

\[
\boxed{
 {\cal U}_5(u,v,x,y,z)=
 \left(u,v,tq,-\frac B2,C\right)
 :\mathbb A^5_{\mathbb Q}\longrightarrow\mathbb A^5_{\mathbb Q}
 }                                                       \tag{3.12}
\]

has

\[
\boxed{\det D{\cal U}_5=1,\qquad
        \operatorname{gdeg}({\cal U}_5)=5.}              \tag{3.13}
\]

Indeed, its Jacobian matrix is block triangular with a `2`-by-`2`
identity block and the vertical determinant-one block of (3.3).  Over the
open `v\pi\ne0`, its complete inverse polynomial at target
`(u,v,\pi,b,c)` is

\[
\boxed{
 E_{u,v,\pi,b,c}(S)
 =v\pi^5S^5+u\pi^4S^4+\pi S^3+bS^2+S-\frac c2.
 }                                                       \tag{3.14}
\]

Let `P(T)` be any monic squarefree quintic over `Q`.  Choose
`a\in\mathbb Q` such that

\[
 P'(a)P^{[3]}(a)\ne0;
\]

only finitely many rational values are excluded.  If

\[
 \frac{P(a+S)}{P'(a)}
 =h_0+S+h_2S^2+h_3S^3+h_4S^4+h_5S^5,
\]

take

\[
\boxed{
 \pi=h_3,\quad b=h_2,\quad c=-2h_0,\quad
 u=\frac{h_4}{h_3^4},\quad v=\frac{h_5}{h_3^5}.
 }                                                       \tag{3.15}
\]

Then (3.14) is exactly `P(a+S)/P'(a)`.  Its discriminant is nonzero,
`v\pi\ne0`, and (3.5) identifies the entire target fiber with

\[
 \operatorname{Spec}\mathbb Q[S]/(P(a+S))
 \simeq \operatorname{Spec}\mathbb Q[T]/(P(T)).          \tag{3.16}
\]

Consequently:

> **Absolute quintic universality.**  One explicit polynomial
> determinant-one map of `A^5_Q`, namely (3.12), has every rank-five
> finite etale `Q`-algebra as a complete target fiber.  In particular, one
> fixed Keller map realizes every quintic number field, so the minimum
> number of fixed Keller maps is one when the ambient dimension is not
> prescribed.

This does not say that the split-seed map of `A^3` realizes every quintic
field.  The two unchanged coordinates in (3.12) are precisely the two
presentation parameters that are fixed in that three-dimensional
specialization.

### 3.3 Exact lifting of parametric quintic families

The construction is compatible with a varying coefficient field, not only
with individual rational polynomials.  Let `B` be an integral rational
variety over `Q`, and let

\[
 P_{\mathbf r}(T)\in\mathbb Q(B)[T]
\]

be a monic separable quintic.  Choose
`a\in\mathbb Q(B)` such that

\[
 P_{\mathbf r}'(a)P_{\mathbf r}^{[3]}(a)\ne0.
\]

After shrinking `B`, all five normalized Hasse coefficients

\[
 \frac{P_{\mathbf r}(a+S)}{P_{\mathbf r}'(a)}
 =h_0+S+h_2S^2+h_3S^3+h_4S^4+h_5S^5
\]

are regular and `h_3h_5` is a unit.  Formula (3.15) defines a rational
target map

\[
 \boxed{
 \kappa_P:B\dashrightarrow\mathbb A^5,\qquad
 \mathbf r\longmapsto
 \left(
 \frac{h_4}{h_3^4},
 \frac{h_5}{h_3^5},
 h_3,h_2,-2h_0
 \right).
 }                                                       \tag{3.17}
\]

Pulling the inverse cover of `{\cal U}_5` back along `\kappa_P` gives

\[
\boxed{
 B\mathop{\times}_{\kappa_P,\mathbb A^5,{\cal U}_5}
 \mathbb A^5
 \simeq
 \operatorname{Spec}_{B}
 \frac{\mathcal O_B[T]}{(P_{\mathbf r}(T))}
 }                                                       \tag{3.18}
\]

on the common squarefree open.  In particular, the degree-five point
extension and its splitting field are unchanged.  If `P` is a generic
polynomial for a transitive group `G<=S_5`, the rational target variety
`\kappa_P(B)` inside this one fixed Keller map parametrizes every
`G`-extension in that degree-five action.

Two small surfaces are completely explicit.  The classical generic
`S_5` polynomial

\[
 P_{r,s}(T)=T^5+rT^3+sT+s
\]

has `a=0` on `rs\ne0`, and (3.17) becomes

\[
\boxed{
 (r,s)\longmapsto
 \left(0,\frac{s^4}{r^5},\frac rs,0,-2\right).
 }                                                       \tag{3.19}
\]

The inverse polynomial (3.14) is exactly `P_{r,s}(S)/s`.
For Brumer's generic `D_5` polynomial

\[
\begin{aligned}
 P_{r,s}(T)={}&T^5+(s-3)T^4+(r-s+3)T^3\\
              &+(s^2-s-2r-1)T^2+rT+s,
\end{aligned}
\]

take `a=0` on `r(r-s+3)\ne0`.  Its target surface is

\[
\boxed{
\begin{aligned}
 u&=\frac{(s-3)r^3}{(r-s+3)^4},&
 v&=\frac{r^4}{(r-s+3)^5},\\
 \pi&=\frac{r-s+3}{r},&
 b&=\frac{s^2-s-2r-1}{r},&
 c&=-\frac{2s}{r}.
\end{aligned}
}                                                        \tag{3.20}
\]

Substitution in (3.14) gives `P_{r,s}(S)/r`.

The remaining three surfaces can also be written exactly.  For any monic
quintic

\[
 P(T)=T^5+a_4T^4+a_3T^3+a_2T^2+a_1T+a_0
\]

with `a_1a_3\ne0`, the origin compiler is

\[
\boxed{
 \kappa(P)=
 \left(
 \frac{a_4a_1^3}{a_3^4},
 \frac{a_1^4}{a_3^5},
 \frac{a_3}{a_1},
 \frac{a_2}{a_1},
 -\frac{2a_0}{a_1}
 \right),
 \qquad
 E_{\kappa(P)}(S)=\frac{P(S)}{a_1}.
}                                                        \tag{3.21}
\]

For Lecacheux's generic `F_{20}` polynomial, put `d=s^2+4` and

\[
\begin{aligned}
 a_4&=td-2s-\frac{17}{4},&
 a_3&=3td+d+\frac{13}{2}s+1,\\
 a_2&=-td-\frac{11}{2}s+8,&
 a_1&=s-6,&a_0&=1.
\end{aligned}
\]

Then (3.21), on `(s-6)a_3\ne0`, is an explicit rational `F_{20}` target
surface, and its inverse polynomial is Lecacheux's polynomial divided by
`s-6`.

For the generic `A_5` polynomial, define

\[
\begin{aligned}
 \Gamma&=5A^2-B^2+3,\\
 \Lambda&=B\Gamma^2-52B\Gamma+576B
          -10\Gamma^2+360\Gamma-3456,\\
 \sigma&=\frac{125\Gamma^2}{4\Lambda},&
 \tau&=\frac{3125\Gamma^5}{256\Lambda^2}.
\end{aligned}
\]

Buhler's polynomial is `P_{A_5}(T)=T^5+\sigma T^3+\tau T+\tau`.
Its Keller target surface simplifies to

\[
\boxed{
 \kappa_{A_5}(A,B)=
 \left(
 0,\frac{3125\Gamma^{10}}{2^{22}\Lambda^3},
 \frac{64\Lambda}{25\Gamma^3},0,-2
 \right),
 \qquad
 E_{\kappa_{A_5}}(S)=\frac{P_{A_5}(S)}{\tau}.
}                                                        \tag{3.22}
\]

Finally, an explicit Hashimoto--Tsunogai `C_5` surface is obtained as
follows.  Put

\[
\begin{aligned}
 Q={}&-A+1+B^2A+7B^2,\\
 R_4={}&A^3+A^2+10B^2A-3A+20B^2+3,\\
 R_3={}&-24B^2A+28B^2+210B^4A+3-28B^2A^2-40B^4-625B^6\\
       &-8A-135B^4A^2-3A^4+2A^5-7B^2A^4+7A^2+44A^3B^2,\\
 R_2={}&4A^4-1+A^6+6A+305B^4+1250B^6+44B^2A^2-220B^4A\\
       &-52A^3B^2+345B^4A^2-2A^5+12B^2A+31B^2A^4
         +11B^2-6A^2,\\
 R_1={}&2A^5-2A^4-8B^2A^4+36A^3B^2-145B^4A^2+3A^2\\
       &-22B^2A^2+4B^2A+120B^4A-2A-13B^2-180B^4-625B^6,\\
 R_0={}&A^3+A^2+7B^2A-B^2.
\end{aligned}
\]

The generic polynomial is

\[
 P_{C_5}(T)
 =T^5-\frac{R_4}{Q}T^4+\frac{R_3}{Q^2}T^3
   +\frac{R_2}{Q^2}T^2+\frac{R_1}{Q^2}T-\frac{R_0}{Q},
\]

and (3.21) gives, on `QR_1R_3\ne0`,

\[
\boxed{
 \kappa_{C_5}(A,B)=
 \left(
 -\frac{QR_4R_1^3}{R_3^4},
 \frac{Q^2R_1^4}{R_3^5},
 \frac{R_3}{R_1},
 \frac{R_2}{R_1},
 \frac{2QR_0}{R_1}
 \right).
}                                                        \tag{3.23}
\]

Its inverse polynomial is `Q^2P_{C_5}(S)/R_1`.  The generic-polynomial
claims and formulas are recorded in
[Jensen--Ledet--Yui, Theorems 2.3.6--2.3.7](https://library.slmath.org/books/Book45/files/book45.pdf),
[Hashimoto--Tsunogai](https://doi.org/10.3792/pjaa.79.142), and
[Kida--Renault--Yokoyama](https://doi.org/10.1142/S1793042109002250).

Thus all five transitive quintic groups now have displayed rational target
surfaces inside the same fixed Keller map.  This is a multi-parameter
fiber-parametric statement.  It is not a one-parameter `G`-parametric
extension, and the ambient map `{\cal U}_5` still has generic monodromy
`S_5`.

## 4. The stack diagram

Let `\mathscr I_N` denote the Keller-fiber incidence stack whose objects
over a scheme `S` are:

- a relative polynomial map `F:A^3_S -> A^3_S` with vertical Jacobian one;
- a target section `y:S -> A^3_S`;
- the condition that `F^{-1}(y)` is finite etale of rank `N`.

One may retain coordinates, or quotient by a specified class of
source--target equivalences; the choice changes the descent question below.
There is always a forgetful morphism

\[
 \rho:\mathscr I_N\longrightarrow BS_N,\qquad
 (F,y)\longmapsto F^{-1}(y).                            \tag{4.1}
\]

The construction gives a commutative triangle

\[
\begin{array}{ccc}
 I_N&\xrightarrow{\ \kappa_N\ }&\mathscr I_N\\
 &\searrow q_N&\downarrow\rho\\
 &&BS_N ,
\end{array}                                             \tag{4.2}
\]

where `q_N` is classified by the finite-etale algebra in (3.5).
The map `q_N` is an etale-local presentation atlas: after a rank-`N`
finite-etale algebra is split, choose `N` distinct rational constants as a
primitive coordinate and choose an origin avoiding the two jet divisors.

The ordered-collision construction itself is choice-free.  For every
rank-`N` finite-etale cover `E->S`, put

\[
 \operatorname{Off}_2(E/S)
 =(E\times_SE)\setminus\Delta_E.
\]

This is a rank-`N(N-1)` finite-etale cover and defines the stack morphism

\[
\boxed{
 \operatorname{Off}_2:BS_N\longrightarrow BS_{N(N-1)}}              \tag{4.3}
\]

induced by the action of `S_N` on ordered distinct pairs.  Pulling (4.3)
back along `q_N` gives exactly `Spec(B_N^{\mathrm{off}})` from (1.14).
Thus the collision cover descends through abstract finite-etale moduli even
though the surrounding quadratic-gauge Keller presentation need not.

Diagram (4.2) is the correct immediate output.  It does **not** reverse to
a canonical arrow

\[
 BS_N\longrightarrow\mathscr I_N.                      \tag{4.4}
\]

Such an arrow would be a choice-free lift of a finite-etale algebra to a
Keller incidence object.  Descent of the displayed construction would
require an isomorphism between the two pullbacks of `\kappa_N` to

\[
 I_N\times_{BS_N}I_N                                   \tag{4.5}
\]

satisfying the cocycle condition.

For `N>=5`, the quadratic-gauge family fails this test generically after
passing even to stable polynomial left--right classes.  The canonical
[generic Tschirnhaus non-descent theorem](GENERIC_TSCHIRNHAUS_NON_DESCENT.md)
works on the clean primitive-presentation groupoid.  The equal-boundary
relation has codimension `N-4`, while the projective Tschirnhaus locus has
codimension `N-3`; outside their union the two presentations define the same
finite-etale algebra but different stable Keller maps.  Hence:

\[
\boxed{\text{the quadratic-gauge atlas map does not descend through full
Tschirnhaus equivalence for }N\ge5.}                    \tag{4.6}
\]

This rules out descent of this particular atlas construction.  It does not
rule out a different, genuinely presentation-free morphism (4.4).

The [rank-three collision-framed descent audit](RANK_THREE_COLLISION_DESCENT.md)
now resolves the finite inertia and presentation cocycle in degree three.
There the off-diagonal pair cover is the full `S_3` frame torsor, and the
unique projective transition between framed root triples lifts to the
foundational factorization map after target localization.  The remaining
rank-three question is global polynomial extension across the explicit
normalizing-denominator divisor; only the scaling torus is denominator-free
within that canonical transport.  The
[rank-four collision-frame audit](RANK_FOUR_COLLISION_CROSS_RATIO.md)
identifies `Conf_3` with the full `S_4` frame torsor and computes the exact
cross-ratio defect
`q_2^2-q_1q_3+q_2q_3e_1+q_3^2e_2`.  Canonical projective root transport
exists only on its zero locus.  The
[rank-four nonprojective lift theorem](RANK_FOUR_NONPROJECTIVE_KELLER_LIFT.md)
now treats the first transverse direction formally: it separates the
ground-field fifth-power twist, gives a rational arithmetic-neutral witness,
reduces the problem to two fibers of one fixed quartic map, and constructs
the exact first-order and all-finite-order formal marked lift.  Its straight
fixed-map target line is finite etale at both endpoints but has the wrong
collision-frame sheet partition.  The `-4` iterate has only two affine
fiber points, so fiber-orbit invariance rules out a polynomial lift of the
straight target translation.  Global polynomial transport remains open only
from target degree nineteen onward with the prescribed frame permutation:
the prime discriminant has degree thirteen, every lower-degree target
self-equivalence is in `mu_5`, and its endpoint orbit fails.  Exact
logarithmic-boundary ranks and a Singular unit-ideal certificate exclude
endpoint degrees thirteen through eighteen, so the open range starts at
degree nineteen.  Stable coarse classes alone do not settle either global
question.

The
[all-rank collision-projective theorem](ALL_RANK_COLLISION_PROJECTIVE_DESCENT.md)
gives the uniform boundary between these cases.  If `r` and `q(r)` are two
primitive coordinates in the normalized degree-`N` quotient, projective
transport exists exactly when the columns `1,r,q(r),r*q(r)` have rank at
most three.  In the basis `1,r,...,r^(N-1)`, after multiplying the last
column by the leading coefficient, these columns form the polynomial matrix

\[
 C_N=
 \begin{pmatrix}
 1&0&q_0&-a_0q_{N-1}\\
 0&1&q_1&a_Nq_0-q_{N-1}\\
 0&0&q_2&a_Nq_1-a_2q_{N-1}\\
 \vdots&\vdots&\vdots&\vdots\\
 0&0&q_{N-1}&a_Nq_{N-2}-a_{N-1}q_{N-1}
 \end{pmatrix}.                                      \tag{4.7}
\]

Here `(a_0,1,a_2,...,a_N)` are the coefficients of the normalized inverse
equation.  Substitution of (3.6) puts every `4`-by-`4` minor of (4.7)
directly in the actual `({\bf u},pi,b,c)` coordinates.  On the primitive
overlap these minors cut out a smooth codimension-`N-3` locus.  Therefore
`Conf_(N-1)` removes all finite `S_N` inertia, but a full presentation
descent must still lift `N-3` independent nonprojective Tschirnhaus
directions away from this locus.

### 4.1 Quantitative seed-descent defect

The projective residual count and the stable Keller-boundary count are
different invariants.  On the compiler slice `a_3=1`, the `N-3` unchanged
seed coordinates are

\[
 (u_4,\ldots,u_N)\in\mathbb G_m^{N-3}.
\]

The residual source--target scaling from the
[quadratic-gauge stable-moduli theorem](QUADRATIC_GAUGE_STABLE_MODULI.md)
acts by

\[
 u_j\longmapsto\alpha^{j+1}u_j.                       \tag{4.8}
\]

After this one redundancy is removed, the intrinsic reconstruction boundary
retains the saturated quotient coordinates

\[
\boxed{
 \Phi_N=
 \left(
 \frac{u_5^5}{u_4^6},
 \frac{u_4u_6}{u_5^2},
 \frac{u_5u_7}{u_6^2},
 \ldots,
 \frac{u_{N-2}u_N}{u_{N-1}^2}
 \right)\in\mathbb G_m^{N-4}.
}                                                       \tag{4.9}
\]

These characters generate the invariant Laurent ring.  By stable
normalization and Fitting functoriality, two compiler points with different
values of `Phi_N` cannot define stably polynomially left--right equivalent
quadratic-gauge maps.  Therefore:

\[
\boxed{
\begin{array}{l}
\text{the framed root embedding has }N-3\text{ independent projective
residuals,}\\
\text{while the quadratic-gauge Keller lift has }N-4\text{ independent
stable boundary moduli.}
\end{array}}                                            \tag{4.10}
\]

In particular, after the single residual scaling, every remaining seed
direction is visible on the intrinsic boundary.  This is the quantitative
form of the obstruction in (4.6).  It still concerns this atlas: it does not
exclude another Keller lift whose presentation cocycle uses different
boundary data.

The two codimensions in (4.10) must not be conflated.  Equality of all
coordinates of `Phi_N` is a codimension-`N-4` condition on a pair of
presentations, whereas projectivity is codimension `N-3`.  Thus equality of
the stable boundary fingerprint is not a projectivity criterion: the generic
equal-fingerprint pair is still nonprojective.  The generic non-descent
statement uses the complement of both loci.

The exact ranks-five-through-seven regression takes the split root
coordinate `r=(1,...,N)` and the primitive Tschirnhaus coordinate
`q(r)=r+r^2`.  The quotient algebras are explicitly isomorphic, but the
applicable coordinates of `Phi_N` all change.  Thus the fibre correspondence
descends while the surrounding stable Keller incidence does not.  The
[bridge card](KELLER_TSCHIRNHAUS_DESCENT_567.md) records the arithmetic
complexity comparison, and its checker pins the exact rational values.

## 5. Braid and symmetric monodromy

Over `C`, the space of monic squarefree degree-`N` polynomials is the
unordered configuration space of `N` points in the affine line.  Its
topological fundamental group is the braid group `B_N`, and the universal
root cover has the standard permutation representation

\[
 B_N\twoheadrightarrow S_N.                            \tag{5.1}
\]

The finite-etale cover in (1.9), or equivalently (3.5), is the pullback of
that universal root cover.  Therefore the monodromy of the marked Keller
fiber is exactly root monodromy; reconstruction adds no new permutation.
Algebraically, over the function field of `U_N` the generic polynomial
still has Galois group `S_N`, because adjoining `a` and localizing at the
discriminant and jet factors do not change the generic splitting group.
Thus

\[
 \boxed{\operatorname{Mon}(q_N)=S_N.}                  \tag{5.2}
\]

On the collision cover this monodromy acts diagonally on ordered pairs:

\[
 [N]\times[N]
 =\{(i,i):i\in[N]\}\sqcup\{(i,j):i\ne j\}.             \tag{5.3}
\]

These are exactly the two factors in (1.14).  The off-diagonal action is
transitive, and the stabilizer of `(1,2)` is the pointwise stabilizer
`S_{N-2}` of the remaining labels.  If `K_N=Frac(R_N)` and `M_N` is the
splitting field of the generic polynomial, then

\[
\boxed{
 \operatorname{Frac}(B_N^{\mathrm{off}})
 =K_N(\alpha_1,\alpha_2)
 =M_N^{S_{N-2}},\qquad
 [K_N(\alpha_1,\alpha_2):K_N]=N(N-1).}                \tag{5.4}
\]

For `N=3`, the stabilizer `S_1` is trivial, so the generic off-diagonal
sheet is the full splitting field:

\[
 B_3^{\mathrm{off}}\otimes_{R_3}K_3\simeq M_3,\qquad
 \operatorname{Gal}(M_3/K_3)\simeq S_3.               \tag{5.5}
\]

This is precisely why a separable nonnormal cubic collision sheet is its
`S_3` normal closure.  For every `N>=4`, `S_{N-2}` is nontrivial and
nonnormal in `S_N`; the generic pair sheet is connected but not Galois and
must not be called the normal closure.  Specialization can lower the
monodromy or disconnect the off-diagonal algebra.

There is a canonical higher extension.  For `1<=m<=N`, let

\[
 \operatorname{Conf}_m(X_N/U_N)
 =X_N^m\setminus\bigcup_{i<j}\{T_i=T_j\},             \tag{5.6}
\]

where `X_N=Spec(B_N)`.  It is the clopen finite-etale cover of ordered
`m`-tuples of distinct roots.  The `S_N`-action is transitive with stabilizer
`S_{N-m}`, and hence

\[
\boxed{
 \operatorname{rk}\operatorname{Conf}_m
 =\frac{N!}{(N-m)!},\qquad
 K(\operatorname{Conf}_m)=M_N^{S_{N-m}}.}             \tag{5.7}
\]

At `m=N-1` and `m=N` the stabilizer is trivial, so the tower reaches the
full `S_N` splitting field in every rank.  The cubic pair sheet is the first
case because `m=2=N-1` when `N=3`.

The extra jet-complement meridians enlarge the fundamental group of the
chosen presentation open, but their permutation action factors through the
usual braid action.  Determining the kernel, including its interaction with
the quadratic-gauge boundary monodromy, is a separate topological problem.

## 6. Essential-dimension lower bound

An atlas dimension is not the relevant invariant for universal
classification: `Spec(Q)->BS_N` is already an etale atlas.  The relevant
condition is *versality over fields*.

Let `X` be an integral finite-type `k`-scheme and suppose a Keller incidence
family over `X` has a rank-`N` finite-etale fiber whose classifying morphism

\[
 X\longrightarrow BS_N
\]

is versal.  By the definition of essential dimension,

\[
 \boxed{\dim X\ge\operatorname{ed}_k(S_N).}             \tag{6.1}
\]

This lower bound applies to every universal Keller parameter space whose
fiber family is versal, independently of the quadratic-gauge formulas.
Conversely, compressing the `S_N`-torsor alone does not automatically
compress the Keller map: a Tschirnhaus compression must also carry descent
data for the surrounding map and target.

For context, Edens--Reichstein identify `ed_k(S_N)` with the minimum
parameter count for the general degree-`N` polynomial under Tschirnhaus
transformations and, in characteristic zero, give

\[
 \left\lfloor\frac{N+1}{2}\right\rfloor
 \le \operatorname{ed}_k(S_N)\le N-3
 \qquad(N\ge6).                                        \tag{6.2}
\]

The exact value is not known in general for `N>=8`; see
[Essential dimension of symmetric groups in prime characteristic,
arXiv:2308.10096](https://arxiv.org/abs/2308.10096).

Equations (3.11) and (6.1) answer different minimization problems:

- `N-3` is achieved and sharp for one relative quadratic-gauge family
  dominant over all normalized polynomial presentations;
- `ed_k(S_N)` is the unavoidable lower bound for a field-versal abstract
  finite-etale family;
- descent from the first object to a parameter space approaching the second
  is obstructed by presentation-dependent stable Keller invariants for
  `N>=5`.

The
[ranks-five-through-seven bridge card](KELLER_TSCHIRNHAUS_DESCENT_567.md)
uses this distinction to propose two further minimization problems: Keller
target dimension and fixed-ambient Keller coordinate degree.  At present
they are organizing definitions, not computed invariants.  The promoted map
only gives the upper bound `ktdim_k(S_N)<=N`, while (6.1) gives the lower
bound `ed_k(S_N)<=ktdim_k(S_N)`.

## 7. What is proved and what remains

Proved here:

1. the literal relative map (1.5), its block-triangular Jacobian, and its
   universal finite-etale marked fiber;
2. the universal ordered collision algebra, its explicit diagonal idempotent,
   its rank-`N(N-1)` off-diagonal factor, and the higher distinct-root tower;
3. the exact compression `U_N=V_N x A^1`;
4. the sharp `N-3`-parameter relative map (3.3);
5. one absolute determinant-one map of `A^N` in every rank `N>=3`, universal
   for rank-`N` finite-etale fibers and absolutely and stably atomic;
6. the stack triangle (4.2), canonical ordered-pair descent, symmetric
   monodromy, and the
   essential-dimension lower bound for versal parameter schemes;
7. generic failure, for `N>=5`, of descent of the quadratic-gauge atlas
   through full presentation equivalence in stable Keller moduli, together
   with the quantitative `N-4`-coordinate boundary fingerprint (4.9), as
   packaged in the canonical generic non-descent theorem;
8. the full `S_N` frame `Conf_(N-1)` and the intrinsic determinantal
   projective locus of smooth codimension `N-3`, including its polynomial
   equations in the universal Keller chart.

Still open:

1. whether a different enhancement of Keller incidence admits a
   presentation-free section of (4.1), or whether the intrinsic Fitting
   decoration can be incorporated into a finite-type enhanced receiver;
2. global polynomial extension of the target-localized rank-three
   projective transport, an all-rank lift on the determinantal projective
   locus, and a rank-four target symmetry of degree at least nineteen
   realizing the formally lifted nonprojective direction with its
   prescribed collision-frame permutation;
3. the kernel of the braid action after the two jet divisors are removed;
4. whether an `ed_k(S_N)`-dimensional versal finite-etale family admits any
   compatible relative polynomial Keller lift;
5. global minimality of the proposed `ktdim` and `kdeg` complexity
   invariants, for which the `N-4` vertical defect of this compiler gives no
   lower bound.

## 8. Exact regression

Run

```bash
.venv/bin/python scripts/verify_universal_relative_keller_map.py
```

The checker verifies the vertical determinant-one identities in the compact
chart, the unchanged-coordinate block promotion in every tested rank, the
normalized inverse equation and coefficientwise compiler through degree
twelve, the `U_N=V_N x A^1` reconstruction, and the compressed
specialization formulas.  It also verifies the adversarial witness cards and
the all-rank `T^N-T-1` target formula through rank twelve, as well as the
degree-drop, bad-translation, and repeated-root boundaries.  On the Osada
`S_N` witnesses in ranks three through eight, it additionally verifies the
two-variable collision presentation, divided-difference identity, Bezout
idempotent, diagonal/off-diagonal Chinese-remainder factors, and exact ranks
`N^2`, `N`, and `N(N-1)`.  An independent permutation audit checks every
ordered distinct `m`-tuple orbit and stabilizer for `1<=m<=N<=8`.

The Lean module `CollisionFiber.lean` formalizes the presentation-independent
tensor collision algebra, diagonal multiplication, obstruction kernel,
ordered-pair functor, and fieldwise obstruction rank.  The Lean
modules `UniversalPromotedBlock.lean`, `UniversalPromotedMap.lean`,
`UniversalPromotedGauge.lean`, `UniversalParameterCompiler.lean`,
`UniversalParameterQuotient.lean`, and
`UniversalParameterWitnesses.lean` formalize the abstract block determinant,
the literal promoted map on an `N`-element coordinate type and its actual
Jacobian-one identity, the coefficient compiler, its quotient-algebra
realization, and three quartic cards.  Symmetric monodromy and
primitive-monodromy atomicity are not yet formalized in Lean and are not
conclusions of a bounded symbolic computation.  The explicit relative
divided-difference splitting and the higher configuration tower are also not
yet Lean theorems.

The separate command

```bash
.venv/bin/python scripts/verify_rank_three_collision_descent.py
```

checks the exact rank-three frame torsor, projective interpolation and
cocycle, quadratic Tschirnhaus boundary ledger, target-localized
factorization transport, and global scaling-torus endpoint.  It does not
classify nonlinear polynomial self-equivalences outside that transport.

The rank-four continuation is checked by

```bash
.venv/bin/python scripts/verify_rank_four_collision_cross_ratio.py
```

It verifies that ordered triples give the full `S_4` frame, factors both the
fourth-root interpolation residual and the cross-ratio difference by the
same explicit defect, separates that hypersurface from the primitive-element
boundary, and writes it in the actual universal-quartic target coordinates.
It does not assert that every Keller equivalence is projective on the root
line.

The all-rank determinantal continuation is checked by

```bash
.venv/bin/python scripts/verify_all_rank_collision_projective_descent.py
```

It verifies the full-frame completion through rank eight, the exact
coefficient matrix (4.7), the cubic and quartic specializations, uniform
projective and nonprojective witnesses, and the `N-3` independent framed
residuals through rank ten.  The uniform theorem is the written
linear-algebra argument; bounded symbolic replay is not substituted for that
proof.

The rank-four nonprojective continuation is checked by

```bash
.venv/bin/python scripts/verify_rank_four_nonprojective_keller_lift.py
```

It verifies the ground-field Kummer class, an arithmetic-neutral primitive
witness, its exact reduction to one fixed quartic map, the finite-etale
straight target line and its sheet mismatch, and both polynomial
first-order lifts.  It also uses the exact two-point fiber at parameter
`-4` to rule out a polynomial lift of the straight target translation.  The
prime discriminant has degree thirteen, so every lower-degree target
self-equivalence is exactly in `mu_5`; its endpoint orbit is also excluded.
The specialized Singular continuation excludes endpoint target degrees
thirteen through eighteen.  The all-finite-order conclusion uses the
separate formal-orbit theorem; no degree-at-least-nineteen endpoint
automorphism is asserted.
