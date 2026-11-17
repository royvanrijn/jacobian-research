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

The collision-algebra and divided-difference interface is credited to
Chloe van der Vlugt's *Collision Ideals and Off-Diagonal Sheets*. The
universal ordered-root application below is the deduction used here.

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
the [finite-etale fibre manuscript](../../../unpublished/common-arithmetic-fibers/main.tex) excludes
full rank two in characteristic zero.  It does not claim that ambient
dimension `N` is minimal.  It also makes no canonical choice of the target
attached to an abstract algebra: the primitive element and the translation
`a` remain presentation choices, exactly as in the atlas discussion of
Section 4.

The hypotheses and conclusions above are those of the displayed
construction. A selected primitive element and translation are presentation
choices, and no canonical choice is supplied by this theorem.

## Current scope and further structure

The relative map, complete fibre identification and absolute atomic map are
the results retained here. The later
[canonical trace-filter receiver](../../jacobian/universal-keller-canonical-filter/PROOF.md)
supplies a different, non-atomic functorial construction.
[Polynomial-profile structure](../../jacobian/universal-keller-structure/PROOF.md)
records the current monodromy and atomicity extension. The published-fibre
and presentation-free descent questions are distinct.

The displayed algebraic arguments are written proofs; the repository curation
does not add formal verification, independent replay or external review.
