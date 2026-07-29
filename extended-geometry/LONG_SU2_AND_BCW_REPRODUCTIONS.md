# Complete local reproductions: `SU(2)`, `SO(3)`, the BCW 79-variable route, and a 21-variable optimization

This note extends the first external-consequences audit with four distinct
outcomes.

1. The `SU(2)` integration formula and Long's lifted witness are now proved
   locally from normalized surface measure on the unit three-sphere.
2. Long's announced `SO(3)` witness is proved locally by pushing Haar measure
   to the third-column sphere.  It is exactly the spherical normalization of
   the repository's two-pair Hopf seed.
3. The complete `3 -> 39 -> 79` Bass--Connell--Wright route is constructed and
   checked exactly.  The fixed-dimensional implication from that map to
   `not GMC(158)` is proved locally in a companion note, following and
   crediting Derksen--van den Essen--Zhao and Zhao.
4. A repository-derived common-factor optimization replaces the conservative
   degree-lowering stage by `3 -> 16`.  The cubic component vector has exact
   rational rank seven, so rank-compressed homogenization needs only 24
   variables.  A later essential-dimension search finds a different
   17-dimensional trace of cubic-output rank six; its 24-dimensional
   homogenization has a three-dimensional constant kernel and gives a
   21-variable cubic-homogeneous collision.  The route-based consequence is
   `not GMC(42)`.
   This is not a formula or dimension claim attributed to Long.

Neither local proof changes the provenance of Christopher D. Long's external
results or constitutes external review.

## 1. A self-contained `SU(2)` integration proof

### 1.1 `SU(2)` as the unit three-sphere

Every element of `SU(2)` has a unique presentation

\[
 g(\alpha,\beta)=
 \begin{pmatrix}
  \alpha&-\overline\beta\\
  \beta&\overline\alpha
 \end{pmatrix},
 \qquad |\alpha|^2+|\beta|^2=1.
\]

Thus `SU(2)` is the unit sphere `S^3` in `C^2=R^4`.  Under the standard
identification with the unit quaternions, left multiplication by
`q=(a,b,c,d)` is represented on `R^4` by

\[
 L_q=
 \begin{pmatrix}
 a&-b&-c&-d\\
 b&a&-d&c\\
 c&d&a&-b\\
 d&-c&b&a
 \end{pmatrix}.
\]

Direct multiplication gives

\[
 L_q^T L_q=(a^2+b^2+c^2+d^2)I_4.
\]

For a unit quaternion, left multiplication is therefore orthogonal.  The
normalized surface measure on `S^3` is left invariant, has total mass one,
and hence is the normalized Haar measure on `SU(2)`.  This uses only the
uniqueness of normalized Haar probability on a compact group.

### 1.2 Exact Hopf-coordinate density

Away from the two measure-zero coordinate circles, put

\[
 \alpha=\sqrt{1-x}\,e^{i\theta_2},
 \qquad
 \beta=\sqrt{x}\,e^{i\theta_1},
 \qquad 0<x<1.
\]

As a map to `R^4`, this is

\[
 \Phi(x,\theta_1,\theta_2)=
 \left(
 \sqrt{1-x}\cos\theta_2,
 \sqrt{1-x}\sin\theta_2,
 \sqrt{x}\cos\theta_1,
 \sqrt{x}\sin\theta_1
 \right).
\]

The three coordinate derivatives are orthogonal and their squared norms are

\[
 \frac{1}{4x(1-x)},\qquad x,\qquad 1-x.
\]

Consequently the Gram determinant is `1/4`, so the surface-volume element is

\[
 \frac12\,dx\,d\theta_1\,d\theta_2.
\]

Integrating this density over `0<=x<=1` and both phases gives total area
`2 pi^2`.  Dividing by that total, normalized Haar measure is exactly

\[
 dx\,\frac{d\theta_1}{2\pi}\,\frac{d\theta_2}{2\pi}.
\]

In particular, `x` is uniform on `[0,1]` and the two phases are independent
normalized circle variables.  The orthogonality and Gram calculations are
checked symbolically by
[`verify_long_su2_haar.py`](../scripts/verify_long_su2_haar.py).

### 1.3 Monomial integration

Use Long's coordinate order

\[
 g=\begin{pmatrix}a&c\\b&d\end{pmatrix},
 \qquad
 a=\alpha,\quad b=\beta,\quad
 c=-\overline\beta,\quad d=\overline\alpha.
\]

For nonnegative integers `r,s,t,u`, Hopf coordinates give

\[
 a^r b^s c^t d^u
 =(-1)^t(1-x)^{(r+u)/2}x^{(s+t)/2}
   e^{i(s-t)\theta_1}e^{i(r-u)\theta_2}.
\]

The two circle integrals vanish unless `r=u` and `s=t`.  In the surviving
case, the beta integral gives

\[
 \int_{SU(2)}a^r b^s c^t d^u\,dg
 =(-1)^s\delta_{r,u}\delta_{s,t}
   \int_0^1(1-x)^r x^s\,dx
 =(-1)^s\delta_{r,u}\delta_{s,t}
   \frac{r!s!}{(r+s+1)!}.                         \tag{1.1}
\]

Now define the polynomial torus substitution

\[
 \beta(z_1,z_2,x)=
 ((1-x)z_2,xz_1,-z_1^{-1},z_2^{-1}).              \tag{1.2}
\]

The monomial `a^r b^s c^t d^u` maps to

\[
 (-1)^t(1-x)^r x^s z_1^{s-t}z_2^{r-u}.
\]

Its torus/Beta integral is again exactly (1.1).  By linearity, for every
polynomial `P(a,b,c,d)`,

\[
 \int_{SU(2)}P\,dg
 =\int_0^1\int_{T^2}P(\beta(z_1,z_2,x))
   \frac{dz_1}{2\pi i z_1}\frac{dz_2}{2\pi i z_2}\,dx.             \tag{1.3}
\]

This is a complete local proof of the Müger--Tuset formula used by Long, not
merely a check conditional on that formula.  Müger and Tuset remain the cited
external source for the formula.

### 1.4 Long's witness

For

\[
 F=(1+c)(ad+b),\qquad G=-c,
\]

the substitution (1.2) gives

\[
 F\longmapsto(1-z_1^{-1})((1-x)+xz_1),
 \qquad
 G\longmapsto z_1^{-1}.
\]

The general beta/binomial computation already proved in the
[external-consequences note](EXTERNAL_CONSEQUENCES_AND_PROVENANCE.md) now
applies through the locally proved formula (1.3), yielding for every `n>=1`

\[
 \int_{SU(2)}F^n\,dg=0,
 \qquad
 \int_{SU(2)}GF^n\,dg=\frac{(-1)^{n-1}}{n+1}\ne0.
\]

Thus the displayed `SU(2)` identities are independently reproduced in full.

### 1.5 Long's `SO(3)` witness

On 28 July 2026, Long announced the following pair in the matrix coordinates
\(R=(r_{ij})\) of \(SO(3)\):

\[
 U=r_{13}+ir_{23},\qquad
 V=r_{13}-ir_{23},\qquad
 T=r_{33},
\]
\[
 P=(1+U)\bigl(V-(2+U)T^2\bigr),\qquad Q=U.
\tag{1.4}
\]

The announcement supplied the identities
\[
 \int_{SO(3)}P^m\,dR=0,\qquad
 \int_{SO(3)}QP^m\,dR
 =\frac{4^m(m!)^2}{(2m+1)!}\ne0
 \quad(m\ge1).
\tag{1.5}
\]
At the time of this audit, no paper, arXiv revision, or stable post URL for
the `SO(3)` announcement had been located.  The external provenance is
therefore the author announcement as quoted to the repository on 28 July
2026; the proof below is a local reproduction, not a claim about priority or
external review.

#### Haar measure reduces to the third-column sphere

The map
\[
 \pi:SO(3)\longrightarrow S^2,\qquad R\longmapsto Re_3
\]
sends a rotation to its third column.  The pushforward of normalized Haar
measure is an \(SO(3)\)-invariant probability measure on \(S^2\), hence is
normalized surface measure.  Since \(P,Q\) depend only on that column, their
group integrals are exactly their spherical integrals.

Write a point of \(S^2\) as \((X,Y,T)\) and put
\[
 U=X+iY,\qquad V=X-iY.
\]
Then
\[
 UV+T^2=1.
\tag{1.6}
\]
The height \(T\) is uniform on \([-1,1]\), and conditional on \(T\), the
phase of \(U=\sqrt{1-T^2}e^{i\theta}\) is uniform.  Thus normalized
spherical integration is
\[
 \frac12\int_{-1}^1\frac1{2\pi}\int_0^{2\pi}(\cdots)\,d\theta\,dT.
\tag{1.7}
\]

#### Endpoint contact gives all moments at once

Using (1.6), the displayed polynomial has the Laurent presentation
\[
\begin{aligned}
 P
 &=(1+U)\left(\frac{1-T^2}{U}-(2+U)T^2\right)\\
 &=\frac{1+U}{U}\left(1-T^2(1+U)^2\right).
\end{aligned}
\tag{1.8}
\]
The apparent \(U^{-1}\) is harmless: (1.4) is polynomial, while (1.8) is
used only for phase constant-term extraction.

Define
\[
 H_m(A)=A^m\int_0^1(1-v^2A^2)^m\,dv.
\tag{1.9}
\]
The integrand is even in the height.  Phase extraction in (1.8) gives
\[
 \int_{SO(3)}P^m\,dR=[u^m]H_m(1+u),\qquad
 \int_{SO(3)}UP^m\,dR=[u^{m-1}]H_m(1+u).
\tag{1.10}
\]
After the change of variable \(w=vA\),
\[
 H_m(A)=A^{m-1}J_m(A),\qquad
 J_m(A)=\int_0^A(1-w^2)^m\,dw.
\tag{1.11}
\]
Now \(J_m'(A)=(1-A^2)^m\) has a zero of order \(m\) at \(A=1\).
Consequently, through Taylor degree \(m\) at \(A=1\),
\[
 H_m(A)=J_m(1)A^{m-1}+O((A-1)^{m+1}).
\tag{1.12}
\]
The coefficient of \(u^m\) is therefore zero, while the coefficient of
\(u^{m-1}\) is \(J_m(1)\).  Finally,
\[
\begin{aligned}
 J_m(1)
 &=\int_0^1(1-w^2)^m\,dw\\
 &=\frac{2^m m!}{(2m+1)!!}
 =\frac{4^m(m!)^2}{(2m+1)!}.
\end{aligned}
\tag{1.13}
\]
Equations (1.10)--(1.13) prove both identities in (1.5) for every
\(m\ge1\).

The functions are of finite type in Mathieu's sense.  Each \(r_{ij}\) is a
matrix coefficient of the defining representation, and polynomials in
matrix coefficients lie in the algebra of representative functions, which
is stable under translations with finite-dimensional orbit span.

#### Relation to the repository seed and consequences

Set
\[
 x=U,\qquad y=\frac V2,\qquad t=T.
\]
Then \(t^2+2xy=1\), and the angular polynomial in the
[two-pair Image--Mathieu counterexample](TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md)
is
\[
 p=(1+x)\left(y-\frac12(2+x)t^2\right)=\frac P2,\qquad Q=x.
\tag{1.14}
\]
Thus Long's `SO(3)` formula is not merely analogous to that construction:
it is exactly the same Hopf-sphere identity, with the factor \(2^m\) in the
mixed moment converted by
\[
 \frac{2^m m!}{(2m+1)!!}
 =\frac{4^m(m!)^2}{(2m+1)!}.
\]
Expanded on the sphere,
\[
 P=V+UV-2T^2-3UT^2-U^2T^2,
\]
so \(P/2\) is also Long's five-term Hopf polynomial recorded in the
[Hopf-lift classification](HOPF_LIFT_CLASSIFICATION.md).

This has four precise consequences.

1. The Mathieu conjecture fails already in the right-\(SO(2)\)-invariant
   subalgebra of finite-type functions on \(SO(3)\), equivalently in
   representative functions on the homogeneous sphere \(SO(3)/SO(2)\).
2. Pullback along the double cover \(SU(2)\to SO(3)\) gives a
   center-invariant `SU(2)` counterexample.  Long's earlier small `SU(2)`
   pair in Section 1.4 is not invariant under the central involution and
   did not itself descend to \(SO(3)\).
3. The sufficient Abelian reduction proposed for \(SO(N)\) in Kevin
   Zwart's
   [*On the Mathieu Conjecture for \(SU(N)\) and \(SO(N)\)*](https://arxiv.org/abs/2304.02648)
   must fail at \(N=3\); here the obstruction is visible in its spherical
   factor.
4. This does not improve the repository's SIC, GVC, or GMC dimension
   bounds: (1.14) identifies an already proved angular seed in a new compact
   group.  Nor does failure for one compact group imply failure of the
   Jacobian conjecture; Mathieu's implication runs from the universal
   compact-group statement to the Jacobian conjecture, not conversely.

### 1.6 The reusable circuit and the higher-group lifting gate

The point of Long's example is not limited to the single pair \(F,G\).
At the Abelian level, put

\[
 f(U,T)=(1-T^{-1})((1-U)+UT).                    \tag{1.15}
\]

The uniform identity used above is the first member of a more general
weighted finite-difference mechanism.  For a nonzero polynomial
\(w\in\mathbb C[U]\), write

\[
 w(U)=(1-U)^{t-1}h(U),\qquad h(1)\ne0,
\]

and set \(D=\deg h\).  The
[weighted constant-term note](PRIME_SEPARATING_RADIAL_MOMENTS.md#no-polynomial-interval-density-has-a-mathieu-kernel)
proves, for every \(m\ge D+t\),

\[
 \int_0^1\operatorname{CT}_T(f(U,T)^m)w(U)\,dU=0,
\]

while

\[
 \int_0^1\operatorname{CT}_T(T^{-t}f(U,T)^m)w(U)\,dU
 =(-1)^{m+t}\frac{h(1)}{t\binom{m+t}{t}}\ne0.     \tag{1.16}
\]

Thus, with \(N=D+t\), the fixed pair \(P=f^N\), \(Q=T^{-t}\)
violates the Mathieu--Zhao conclusion for every positive power of \(P\).
For \(w=1\), this is exactly Long's sequence
\((-1)^{m-1}/(m+1)\).  The proof is the same circuit calculation at a
higher finite-difference order: after the beta integral, the Bernstein
coefficient is a polynomial in the angular index, so the full alternating
difference vanishes and a boundary term survives after multiplication by
\(T^{-t}\).

This interacts directly with Müger--Tuset's general integration formula.
For a direct product of compact connected simple groups and a torus, their
finite-type reduction supplies integers \(N_K,M_K\), a nonzero monomial

\[
 \delta(x)=c\prod_{j=1}^{N_K}x_j^{q_j},
 \qquad q_j\ \hbox{odd},
\]

and a linear map

\[
 \mathbb C[K]\longrightarrow
 \mathbb C[x_1,\ldots,x_{N_K},
           z_1^{\mathord\pm1},\ldots,z_{M_K}^{\mathord\pm1}],
 \qquad H\longmapsto\widetilde H,
\]

that preserves every joint moment
\(\int_K H^aJ^b\) after weighted product integration.  Arbitrary compact
connected Lie groups are reached from such a product by a finite central
quotient.

At the level of the **full ambient Laurent-polynomial ring**, (1.16)
therefore applies to every Müger--Tuset monomial density as soon as there is
at least one interval variable and one torus variable.  Indeed, choose one
factor \(x_j^{q_j}\), use (1.15) in \((x_j,z_1)\), ignore the other variables,
and take

\[
 P=f^{q_j+1},\qquad Q=z_1^{-1}.                   \tag{1.17}
\]

For every \(\ell\ge1\), the weighted pure integral of \(P^\ell\) is zero
and the weighted integral of \(QP^\ell\) is nonzero.  The unused factors
contribute only nonzero constants.  This is an exact consequence of the
local polynomial-density theorem, not a bounded search.  It shows that the
strong ambient Abelian conjecture cannot be repaired merely by retaining
the monomial weights occurring in the general Lie-group formula.

It does **not** give a counterexample on every compact connected Lie group.
The group-dependent image

\[
 \mathcal A(K)=\{\widetilde H:H\in\mathbb C[K]\}
\]

is generally a proper subspace of the ambient Laurent-polynomial ring.
Moreover, Müger and Tuset prove the moment identities using a linear
replacement map; they explicitly do not claim that this replacement is an
algebra homomorphism.  A higher-group counterexample must therefore pass
the following lifting gate:

> Find \(H,J\in\mathbb C[K]\) whose Müger--Tuset replacements are the
> circuit pair (1.17), or another pair in \(\mathcal A(K)\) with the same
> all-order pure and mixed moments.

For a finite central quotient there is an additional descent gate: the
chosen functions on the product cover must be invariant under the finite
kernel.  Groups with an actual `SU(2)` direct factor inherit Long's witness
immediately by ignoring the other factors.  A root `SU(2)` subgroup inside
a simple higher-rank group is not enough, because a function on that
subgroup need not extend to a group function whose integration-coordinate
image ignores all other root factors.

This isolates a reproducible higher-rank program.

1. Start with `SU(3)` and a fixed reduced word for the longest Weyl element.
   Compute the replacements of low-weight matrix coefficients and
   generalized minors, then test exact membership of the powered circuit
   (1.17) and its multiplier in their span.
2. If literal membership fails, search inside \(\mathcal A(K)\) for the
   weaker invariant condition that the Bernstein coefficient sequence has
   uniformly bounded degree and a nonzero truncated boundary difference.
3. Only after an all-order identity is obtained on the product cover, test
   invariance under the relevant finite central kernel.

The circuit already has two separate algebraic descendants in the active
Image-Mathieu work.  The
[four-term three-pair witness](THREE_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md)
is a one-pair bihomogeneous lift of Long's `SU(2)` seed, while the
[Dvorsky--Long witness](DVORSKY_GVC5_COUNTEREXAMPLE.md) is a five-variable
constant-coefficient GVC and five-pair SIC lift.  The later
[two-pair counterexample](TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md) proves
that the minimum failing SIC pair dimension is two, but it is a genuinely
nonseparable rank-five construction rather than a Müger--Tuset lift.
Likewise, the Gaussian counterexamples involve factorial radial weights
rather than the beta weights in (1.15)--(1.16).  The present observation
therefore sharpens the compact-group lifting program without changing the
repository's GMC, SIC, or BCW dimension bounds.

The status distinction is therefore:

- the weighted ambient failure is an exact corollary of the locally proved
  polynomial-density theorem;
- the `SU(2)` lift is the reproduced theorem of Long;
- membership in \(\mathcal A(K)\) for a simple higher-rank group remains an
  open problem.

## 2. The exact `3 -> 39 -> 79` BCW route

### 2.1 Starting normalization and support

Long uses the determinant-one presentation `L` which satisfies

\[
 L=\operatorname{diag}(1/2,1/2,-1/2)
   \circ F\circ\operatorname{diag}(1,2,2)
\]

for the repository's foundational determinant `-2` map `F`.  Composing the
target with `(u,v,w) -> (-w,v,u)` gives a map `G` with identity linear part:

\[
 \begin{aligned}
 G_1={}&x-3x^2y-x^3z,\\
 G_2={}&y+3xz+24xy^2+12x^2yz+36x^2y^3+12x^3y^2z,\\
 G_3={}&z+8y^2+6xyz+28xy^3+12x^2y^2z
          +24x^2y^4+8x^3y^3z.
 \end{aligned}
\]

Its terms of degrees `4,5,6,7` occur with exact counts

\[
 3,\quad2,\quad2,\quad1.                            \tag{2.1}
\]

The determinant, normalization, support counts, and three-point collision are
all checked exactly by the route script.

### 2.2 One stable degree-lowering step

Suppose one coordinate `f_i` contains a monomial `ab`, with `deg a=p`,
`deg b=q`, `p,q>=2`.  Adjoin variables `Y,Z` and replace the stabilized map by

\[
 \widetilde f_i=f_i-(Y+a)(Z+b),
 \qquad
 \widetilde f_{n+1}=Y+a,
 \qquad
 \widetilde f_{n+2}=Z+b,                            \tag{2.2}
\]

leaving the other coordinates unchanged.  This is an exact polynomial
left--right equivalence.  Indeed, first apply the source automorphism

\[
 (x,Y,Z)\longmapsto(x,Y+a(x),Z+b(x))
\]

to `f x id`, then apply the target automorphism which subtracts the product of
the last two outputs from output `i`.

The term `ab` cancels.  The new nonlinear terms have degrees at most

\[
 p+1,\quad q+1,\quad p,\quad q,\quad2.              \tag{2.3}
\]

The determinant, generic degree, noninvertibility, and collision schemes are
preserved.  On a collision point `x`, adjoining `Y=-a(x), Z=-b(x)` transports
the collision explicitly.

### 2.3 Eighteen balanced steps

Let `c(d)=0` for `d<=3`.  Factoring a degree-`d` monomial into balanced degrees
`p+q=d` and using (2.3) gives

\[
 c(d)\leq1+c(p+1)+c(q+1)+c(p)+c(q).
\]

The choices

\[
 4=2+2,\quad5=2+3,\quad6=3+3,\quad7=3+4
\]

give the upper bounds

\[
 c(4)\leq1,\qquad c(5)\leq2,\qquad
 c(6)\leq3,\qquad c(7)\leq5.
\]

Together with (2.1), the total is

\[
 3c(4)+2c(5)+2c(6)+c(7)\leq18.
\]

The exact implementation performs the degree sequence

\[
 7,6,6,5,5,5,\underbrace{4,\ldots,4}_{12\text{ times}},
\]

and terminates with a noninjective determinant-one map

\[
 K(X)=X+Q(X)+C(X),\qquad X\in\mathbb A^{39},        \tag{2.4}
\]

where `Q` is quadratic homogeneous and `C` is cubic homogeneous.  Since each
step adds two variables, the dimension is exactly

\[
 3+2\cdot18=39.
\]

### 2.4 Nilpotent form in 78 variables

On variables `(X,Y)` in `A^39 x A^39`, define

\[
 U(X,Y)=(X+Q(X)+Y,\;Y-C(X))=(X,Y)+N(X,Y).           \tag{2.5}
\]

This map is stably equivalent to `K`.  Starting with `K x id`, apply

\[
 A(X,Y)=(X,Y-C(X))
\]

on the source and then

\[
 P(R,S)=(R+S,S)
\]

on the target.  The result is exactly (2.5).

It remains to prove that `JN` is nilpotent, not merely that `det(I+JN)=1`.
For an indeterminate `t`, put

\[
 E_t(X)=X+tQ(X)+t^2C(X)=t^{-1}K(tX).
\]

Hence `det DE_t=1`.  More explicitly,

\[
 (X,Y)+tN(X,Y)
 =P_t\circ(E_t\times\operatorname{id})\circ A_t(X,Y),
\]

where

\[
 A_t(X,Y)=(X,Y-tC(X)),\qquad P_t(R,S)=(R+tS,S).
\]

Therefore

\[
 \det(I+tJN)=1                                      \tag{2.6}
\]

as a polynomial identity in `t`.  All nonconstant coefficients of the
characteristic polynomial of `JN` vanish, so `JN` is nilpotent.

### 2.5 Cubic homogenization in 79 variables

Adjoin `T` and homogenize the degree-one, degree-two, and degree-three parts of
`N`:

\[
 H(X,Y,T)=(YT^2+Q(X)T,\;-C(X),\;0).
\]

Every component of `H` is cubic homogeneous.  The map

\[
 V(X,Y,T)=(X,Y,T)+H(X,Y,T)                          \tag{2.7}
\]

lives in

\[
 2\cdot39+1=79
\]

variables.  To verify its Jacobian, write `Z=(X,Y)`.  On `T!=0`,

\[
 H(Z,T)=T^3N(Z/T),
 \qquad
 \partial_ZH=T^2JN(Z/T).
\]

Equation (2.6), with `t=T^2`, gives `det DV=1` on this dense open, hence
everywhere by polynomial identity.

At `T=1`, map (2.7) restricts to `(U,1)`.  If `X_i` are the three transported
collision points of `K`, then

\[
 (X_i,C(X_i),1)
\]

are three distinct points with one common image under `V`.  Thus `V` is an
explicit noninvertible cubic-homogeneous Keller map in 79 variables.

The complete construction and collision are checked by
[`verify_long_bcw_79_route.py`](../scripts/verify_long_bcw_79_route.py).
That script writes the exact
[79-variable sparse artifact](../artifacts/generated-results/long_bcw_79_counterexample.json),
which is replayed without SymPy by
[`audit_long_bcw_79_independent.py`](../scripts/audit_long_bcw_79_independent.py).

### 2.6 The final Gaussian implication

For each fixed `r`, the needed implication is

\[
 \mathrm{GMC}(2r)\Longrightarrow\mathrm{SIC}(r).
\]

and `SIC(r)` forces a cubic-homogeneous Keller map in `r` variables to be
invertible.  The companion
[fixed-dimensional proof](FIXED_GMC_SIC_PROOF.md) derives both statements from
Gaussian contraction, an uncountable-field countable-union lemma, and a
coefficient proof of the formal inversion identity.  Applying its
contrapositive to (2.7), with `r=79`, gives

\[
 \neg\mathrm{GMC}(158).
\]

This is now locally proved, but it is still nonconstructive at the Gaussian
witness step because the passage from pointwise thresholds to one uniform
threshold uses a countable-union argument.  The proof reproduces only the
fixed-dimensional implication needed here, not all results of DVEZ or Zhao.

## 3. Shared-factor and rank-compressed optimization: `3 -> 16 -> 24`

The 79-variable construction above remains the exact reproduction of Long's
conservative route.  It is not dimension-minimal: its 18 steps expose two
fresh factors independently even when a factor has already appeared as an
output coordinate.

### 3.1 Reusable-factor elementary equivalence

Suppose a previous source shear has exposed a polynomial `a(x)` as the output

\[
 A=Y+a(x).                                           \tag{3.1}
\]

If coordinate `i` contains `c a(x)b(x)` and `i` is not the `A` coordinate,
adjoin only one new variable `Z`, apply the source shear

\[
 Z\longmapsto Z+b(x),
\]

and then the elementary target shear

\[
 T_i\longmapsto T_i-cA T_Z.                         \tag{3.2}
\]

After the source shear, `T_Z=Z+b(x)`, so (3.2) cancels `cab`.  Both
automorphisms have determinant one.  If `b` was exposed previously too, no
new variable is needed: subtract `cAB` directly.  If `a=b` was not exposed,
one variable suffices—expose `A=Y+a` and subtract `cA^2`.  These are exact
stable left--right equivalences, not arithmetic-circuit substitutions made
outside the map category.

At a transported collision point, set every new source variable to the
negative value of its exposed factor.  All exposed-factor outputs are then
zero, so each target shear preserves the common image exactly.

### 3.2 Deterministic factor search

For each current top-degree monomial, enumerate every unordered factor split
`ab` with both degrees between two and `d-2`.  Retain a registry of outputs of
the form (3.1), and rank candidates lexicographically by

\[
 (\text{new maximum degree},
   \sum_{\deg m>3}(\deg m-3)^2,\
   \text{number of high terms},\
   \text{new variables},\
   -\text{reusable factors}).                        \tag{3.3}
\]

Freezing the best trace from a deterministic width-24 beam search produces 17
target cancellations with degree
sequence

\[
 7,6,6,5,5,5,5,\underbrace{4,\ldots,4}_{10\text{ times}},
\]

but the corresponding new-variable counts are

\[
 2,0,1,0,1,1,1,0,0,1,1,1,1,1,1,0,1.              \tag{3.4}
\]

Their sum is 13.  The resulting determinant-one map therefore has degree at
most three in

\[
 3+13=16                                             \tag{3.5}
\]

variables, identity linear part, and the transported rational three-point
collision.  Five zero-cost cancellations in (3.4) reuse both exposed factors;
the square cancellation at the second step is the first of them.

This is a certified upper bound, not a minimality theorem.  The greedy score
is deliberately recorded so SAT, MILP, dynamic-programming, or beam searches
can seek a still smaller exposure registry without changing the certificate
format.

### 3.3 Rank-compressed cubic homogenization

Write the 16-variable map as `K=X+Q+C`, with `Q,C` homogeneous of degrees two
and three.  Let `k` be the row rank over `Q` of the coefficient matrix of the
component vector `C`.  Choose independent component polynomials
`c=(c_1,...,c_k)` and the unique constant matrix `B` such that

\[
 C(X)=B c(X).                                       \tag{3.6}
\]

Only `Y in A^k` is needed.  Put

\[
 U(X,Y)=(X+Q(X)+BY,Y-c(X)).                         \tag{3.7}
\]

This is stably left--right equivalent to `K`: precompose `K times id` with
`A(X,Y)=(X,Y-c(X))`, then postcompose with `P(R,S)=(R+BS,S)`.  More strongly,
for

\[
 E_t(X)=X+tQ(X)+t^2C(X)=t^{-1}K(tX),
\]

define `A_t(X,Y)=(X,Y-tc(X))` and
`P_t(R,S)=(R+tBS,S)`.  Direct substitution gives

\[
 \operatorname{id}+tN
 =P_t\circ(E_t\times\operatorname{id})\circ A_t,
 \qquad N=(Q+BY,-c).                                \tag{3.8}
\]

All source and target shears have determinant one, and
`det DE_t=det DK(tX)=1`; hence `det(I+tJN)=1`.  Equivalently, the Schur
complement in the homogenized Jacobian is

\[
 I+tJQ+t^2B Jc=DE_t.
\]

Therefore

\[
 V_{n+k+1}(X,Y,T)
 =(X,Y,T)+(TQ(X)+T^2BY,-c(X),0)                     \tag{3.9}
\]

is cubic homogeneous and Keller.  If `K(p)=q`, then
`V(p,c(p),1)=(q,0,1)`, so every collision transports.

For the frozen 16-variable trace, precisely the components numbered
`0,1,2,3,4,6,8` are nonzero.  Exact rational row reduction shows that all
seven are independent.  Thus `k=7`, not merely `k<=7`, and (3.9) has
`16+7+1=24` variables.  Its transported rational three-point collision makes
it noninvertible.  A separate 17-dimensional trace has cubic-output rank six;
its 24-dimensional homogenization has a three-dimensional constant kernel.
Quotienting those directions preserves the collision and gives a 21-variable
cubic-homogeneous Keller map.  The locally proved
fixed-dimensional implication now gives

\[
 \boxed{\neg\mathrm{GMC}(42)}.                     \tag{3.10}
\]

This improves only the nonexplicit route-based dimension bound.  Long's
direct three-real-Gaussian witness remains far stronger and independently
authored.

The SymPy generator
[`verify_shared_bcw_33_route.py`](../scripts/verify_shared_bcw_33_route.py)
writes the exact
[33-variable sparse artifact](../artifacts/generated-results/shared_bcw_33_counterexample.json).
The dependency-free
[`audit_shared_bcw_33_independent.py`](../scripts/audit_shared_bcw_33_independent.py)
replays all factor exposures and target shears from the original map and
reconstructs the cubic collision without importing the generator.
The general rank factorization is implemented in
[`rank_compressed_bcw_homogenization.py`](../scripts/rank_compressed_bcw_homogenization.py).
The generator
[`verify_rank_compressed_bcw_24_route.py`](../scripts/verify_rank_compressed_bcw_24_route.py)
writes the exact
[24-variable sparse artifact](../artifacts/generated-results/rank_compressed_bcw_24_counterexample.json),
and
[`audit_rank_compressed_bcw_24_independent.py`](../scripts/audit_rank_compressed_bcw_24_independent.py)
independently recomputes the rational rank, replays (3.8) with sparse
polynomials, reconstructs the map, and checks the collision using only the
standard library.

The quotient generator
[`verify_constant_kernel_bcw_22_route.py`](../scripts/verify_constant_kernel_bcw_22_route.py)
writes the
[22-variable sparse artifact](../artifacts/generated-results/constant_kernel_bcw_22_counterexample.json).
The dependency-free
[`audit_constant_kernel_bcw_22_independent.py`](../scripts/audit_constant_kernel_bcw_22_independent.py)
parses the 24-variable source artifact, recomputes both kernel vectors, checks
`BK=0`, `BC=I`, and `H=HCB`, reconstructs the quotient, verifies cubic
homogeneity and the descended collision, and checks the block-triangular
determinant factorization.

The frozen generator
[`verify_essential_bcw_21_route.py`](../scripts/verify_essential_bcw_21_route.py)
writes the
[21-variable sparse artifact](../artifacts/generated-results/essential_bcw_21_counterexample.json).
The dependency-free
[`audit_essential_bcw_21_independent.py`](../scripts/audit_essential_bcw_21_independent.py)
replays the 17 cancellations from the original map, independently recovers
cubic-output rank six, reconstructs the 24-dimensional homogenization,
recomputes its three-dimensional constant kernel, and rebuilds the 21D
quotient, collision, and triangular determinant factorization.

### 3.4 What remains after 42

The bound 42 is still an upper bound, not a minimality theorem.  The executable
[`search_rank_aware_bcw.py`](../scripts/search_rank_aware_bcw.py) reconstructs
the monomial-factor beam search, deduplicates exact polynomial states, and
uses `s+rank(C)` only as its partial-state beam heuristic.  Every completed
trace is now rank-compressed, homogenized, constant-kernel quotiented, and
compared by final essential dimension.  The fuller
[`search_essential_bcw.py`](../scripts/search_essential_bcw.py) also transports
the collision and computes good-prime cyclic invariant-row-module diagnostics.
Historically, both the legacy degree-first and genuinely rank-first orderings, at width 128,
finish with `s+rank(C)=20`; neither finds a 23-variable map.  Across the 232
completed traces retained by the degree-first run, no trace passes eight exact
necessary samples of `det(I+sJQ+tJC)=1`.  These are finite search results, not
lower-bound proofs.

A kernel-aware ordering now reranks a bounded prebeam by the modular
post-kernel essential dimension of each partial quadratic--cubic truncation.
At widths 24 and 64 it completes 44 and 112 terminal traces respectively;
every terminal coefficient matrix has modular rank 21, so none can have
essential dimension below 21 over `Q`.  On the width-24 terminal set, no trace
passes the eight exact two-parameter samples; the closest traces fail four.
Again, this is bounded negative search evidence, not minimality.

The present 16-variable map itself definitely cannot avoid doubling.  At
`X=(1,...,1)` its known scaling family still has

\[
 \det(I+JQ+JC)=\det(I+2JQ+4JC)=1,
\]

but the independent specializations give

\[
 \det(I+JC)=-4160,\qquad \det(I+JQ)=-78.            \tag{3.11}
\]

Thus the desired two-parameter identity fails before homogenization.  The
exact checker
[`verify_two_parameter_bcw_obstruction.py`](../scripts/verify_two_parameter_bcw_obstruction.py)
records this obstruction.  Further improvement now requires a broader
equivalence class—most plausibly polynomial-factor reuse and multi-term
cancellation—or a larger search that escapes the width-128 monomial beam.

## Reproduction

Run

```bash
.venv/bin/python scripts/verify_long_su2_haar.py
.venv/bin/python scripts/verify_long_xz_mathieu.py
.venv/bin/python scripts/verify_beta_radial_mathieu_counterexamples.py
.venv/bin/python scripts/verify_long_bcw_79_route.py
python3 scripts/audit_long_bcw_79_independent.py
.venv/bin/python scripts/verify_shared_bcw_33_route.py
python3 scripts/audit_shared_bcw_33_independent.py
.venv/bin/python scripts/verify_rank_compressed_bcw_24_route.py
python3 scripts/audit_rank_compressed_bcw_24_independent.py
.venv/bin/python scripts/verify_constant_kernel_bcw_22_route.py
python3 scripts/audit_constant_kernel_bcw_22_independent.py
.venv/bin/python scripts/verify_essential_bcw_21_route.py
python3 scripts/audit_essential_bcw_21_independent.py
.venv/bin/python scripts/audit_bcw_21_linear_quotients.py
python3 scripts/audit_bcw_21_affine_vector_symmetries.py
.venv/bin/python scripts/verify_two_parameter_bcw_obstruction.py
python3 scripts/verify_fixed_gmc_sic_bridge.py
```

The first two scripts jointly certify the complete `SU(2)` proof; the second
also directly checks the `SO(3)` spherical moments through order fifteen and
the endpoint-jet identities through order one hundred.  The third checks the
integer-Beta and arbitrary polynomial-density extensions of the same
circuit.  The next two construct and independently replay Long's
conservative 79-variable route.  The next pair record and replay the
repository's shared-factor baseline, and the following pair construct and
independently replay its rank-compressed 24-variable homogenization.  The
next pair construct and independently replay its 22-variable constant-kernel
quotient.  The following pair freeze and independently replay the improved
21-variable essential quotient.  The next two scripts exclude further
collision-preserving linear quotients and affine-vector-field translation
symmetries.  The next script checks the two-parameter shortcut obstruction;
the final script checks the coefficient skeleton of the fixed-dimensional
implication.  None of these replaces the repository's separate 95-variable
cubic-homogeneous artifact.

The external theorem inputs and construction sources are:

- Christopher D. Long, [*Small Counterexamples to the Gaussian Moments
  Conjecture*](https://arxiv.org/abs/2607.18186), arXiv:2607.18186v1;
- Christopher D. Long, [*Counterexamples to the (xz)-Conjecture and the
  Mathieu Conjecture for (SU(2))*](https://arxiv.org/abs/2607.19012),
  arXiv:2607.19012v1;
- Christopher D. Long, public `SO(3)` counterexample announcement,
  28 July 2026, quoted in Section 1.5; no archival URL was located at audit
  time;
- Kevin Zwart, [*On the Mathieu Conjecture for \(SU(N)\) and
  \(SO(N)\)*](https://arxiv.org/abs/2304.02648), arXiv:2304.02648;
- Michael Müger and Lars Tuset, [*An integral formula for Lie groups, and
  the Mathieu conjecture reduced to Abelian non-Lie
  conjectures*](https://arxiv.org/abs/2410.11622), arXiv:2410.11622v2,
  accepted in Advances in Mathematics;
- Hyman Bass, Edwin H. Connell, and David Wright,
  [*The Jacobian Conjecture: Reduction of Degree and Formal Expansion of the
  Inverse*](https://doi.org/10.1090/S0273-0979-1982-15032-7), Bulletin of the
  AMS 7 (1982), 287--330;
- Harm Derksen, Arno van den Essen, and Wenhua Zhao,
  [*The Gaussian Moments Conjecture and the Jacobian
  Conjecture*](https://arxiv.org/abs/1506.05192), Israel Journal of Mathematics
  219 (2017), 917--928; and
- Wenhua Zhao, [*Images of commuting differential operators of order one with
  constant leading coefficients*](https://arxiv.org/abs/0902.0210), Journal
  of Algebra 324 (2010), 231--247.
