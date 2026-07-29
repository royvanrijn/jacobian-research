# Local moduli of the bidegree-\((4,4)\) two-pair counterexample

## 1. Scope and conventions

This note studies only the displayed point

\[
 F=(R+Z)\left(R^2W-\frac12(2R+Z)T^2\right)
 \in V_4=\operatorname{Sym}^4(U^*)\otimes\operatorname{Sym}^4(U),
 \tag{1.1}
\]

where

\[
 R=\xi _1z_1+\xi _2z_2,\quad Z=\xi _1z_2,\quad
 W=2\xi _2z_1,\quad T=\xi _1z_1-\xi _2z_2
 \tag{1.2}
\]

and \(T^2=R^2-2ZW\).  No other rank stratum and no sparse or global
search is used.

The contraction-preserving linear group is

\[
 G=\operatorname{GL}(U)/\mathbb G_m\simeq\operatorname{PGL}_2,
 \tag{1.3}
\]

acting by \(z\mapsto gz\), \(\xi\mapsto g^{-T}\xi\).  The center of
\(\operatorname{GL}(U)\) acts trivially.  Write

\[
 \ell_n(H)=\mathcal E_2(HF^n)\qquad(n\geq0).
 \tag{1.4}
\]

The all-moment-zero scheme is denoted by
\[
 {\cal Z}=V(\mathcal E_2(X),\mathcal E_2(X^2),\ldots)\subset V_4.
 \tag{1.5}
\]

## 2. Stabilizer and orbit

On the determinant-one cover of \(G\), put
\[
 g=\begin{pmatrix}a&b\\c&d\end{pmatrix}.
\]
Exact coefficient comparison in \(gF=F\), together with \(ad-bc=1\),
has reduced Gröbner basis
\[
 b,\quad c,\quad a-d,\quad(d-1)(d+1).
 \tag{2.1}
\]
Thus
\[
 \operatorname{Stab}_{\operatorname{SL}_2}(F)=\{\pm I\},\qquad
 \boxed{\operatorname{Stab}_{G}(F)=1}.
 \tag{2.2}
\]
In particular the contraction-preserving orbit is smooth and
\[
 \boxed{\dim(GF)=3}.
 \tag{2.3}
\]

## 3. The exact all-order tangent space

Use the ordered monomial basis
\[
 e_{ij}=\xi _1^i\xi _2^{4-i}z_1^jz_2^{4-j},
 \qquad 0\leq i,j\leq4.
 \tag{3.1}
\]
Differentiating every moment at \(F\) gives
\[
 d(\mathcal E_2(X^m))_F(H)=m\ell_{m-1}(H).
 \tag{3.2}
\]
The infinite ideal is not replaced by a guessed finite prefix.  Instead,
the same Hopf coefficient extraction used for the counterexample gives,
for every basis monomial,
\[
\begin{split}
 \frac{\ell_n(e_{ij})}{(4n+5)!}
 ={}&2^{-n-j-4+i}
 \sum_{k=0}^{n}(-1)^k\binom nk
 \binom{n+2k}{\,n-i+j\,}\\
 &\mathrel{\phantom{=}}\cdot
 \frac12\int_{-1}^{1}
 t^{2k}(1+t)^j(1-t)^{4-i}\,dt .
\end{split}
\tag{3.3}
\]
This is an identity for every \(n\), not a bounded moment check.

Expanding the beta integral and applying the same change
\(X=1+x,\ v=tX\) as in the seed proof gives the following tail identity.
For \(n\geq16\), as functionals on \(V_4\),
\[
\ell_n\ \sim\
\ell_{12}
 -\frac{n-12}{1521520(n-13)}\ell_{13}
 +\frac{n-12}{9199596806400(n-14)}\ell_{14}
 -\frac{n-12}{216992729791918080000(n-15)}\ell_{15},
\tag{3.4}
\]
where \(\sim\) means equality up to a nonzero scalar depending only on
\(n\).  The scalar obeys
\[
 \frac{s_{n+1}}{s_n}
 =\frac{8(n-12)(n+1)(n+2)(4n+7)(4n+9)}{n-15}.
\tag{3.5}
\]
Direct reduction of (3.3), still symbolically in \(n\), puts
\(\ell_{12},\ldots,\ell_{15}\) in the span of
\(\ell_0,\ldots,\ell_{11}\).  Those first twelve functionals are
independent.  Hence (3.3)--(3.5), rather than a finite-prefix
assumption, prove
\[
\boxed{
 T_F{\cal Z}
 =\bigcap_{n\geq0}\ker\ell_n
 =\bigcap_{n=0}^{11}\ker\ell_n,\qquad
 \dim T_F{\cal Z}=13.}
\tag{3.6}
\]
The exact \(12\) by \(25\) matrix and an exact \(25\) by \(13\) kernel
basis are stored by the checker.

The orbit tangent \(\mathfrak {sl}_2F\) has dimension three and is
contained in (3.6).  The radial direction \(kF\) is a fourth independent
direction in (3.6).

## 4. There is no nullcone tangent at \(F\)

Identify \(V_4\) with \(\operatorname{End}(\operatorname{Sym}^4U)\).
If \(C=(c_{ij})\) is the coefficient matrix in (3.1), put
\[
 D=\operatorname{diag}(0!4!,1!3!,2!2!,3!1!,4!0!),\qquad
 I_2(C)=\operatorname{tr}((DC)^2).
\tag{4.1}
\]
This is a homogeneous quadratic \(G\)-invariant, and exact evaluation
gives
\[
 \boxed{I_2(F)=1152\ne0.}
\tag{4.2}
\]
The nullcone is contained in \(V(I_2)\), so the principal open
\(D(I_2)\) is a Zariski neighbourhood of \(F\) disjoint from the
nullcone.  Therefore the nullcone has no local germ at \(F\).  In
particular, a Zariski tangent space \(T_FN_4\) is not defined because
\(F\notin N_4\); if tangent directions to the local intersection are
interpreted set-theoretically, the answer is the empty set, not a vector
subspace.

### 4.1 The minimum missing invariant

The quadratic separator in (4.1) is not an accidental trace formula.
It is the first member of a complete low-degree calculation.  Put
\[
 A=C^TD,\qquad
 D=\operatorname{diag}(0!4!,1!3!,2!2!,3!1!,4!0!).
 \tag{4.3}
\]
This is the operator matrix attached to the coefficient tensor; transposing
changes neither its trace powers nor its determinant.  Under conjugation,
\[
 V_4=\operatorname{End}(\operatorname{Sym}^4)
 \cong
 \operatorname{Sym}^0\oplus\operatorname{Sym}^2
 \oplus\operatorname{Sym}^4\oplus\operatorname{Sym}^6
 \oplus\operatorname{Sym}^8.
 \tag{4.4}
\]
Let \(P_{2r}\) be the equivariant projector onto
\(\operatorname{Sym}^{2r}\), obtained from the adjoint Casimir, and set
\[
 A_{2r}=P_{2r}(A),\qquad q_{2r}=\operatorname{tr}(A_{2r}^2).
 \tag{4.5}
\]
The Casimir eigenvalues on the five summands are
\[
 0,\ 4,\ 12,\ 24,\ 40,
 \tag{4.6}
\]
so (4.5) is an exact rational construction.  Multiplicity freeness and
self-duality show that
\[
 (R_4)_2
 =\langle q_0,q_2,q_4,q_6,q_8\rangle,\qquad
 q_0=\frac{\mu_1^2}{5}.
 \tag{4.7}
\]
At the displayed \(F\), exact projection gives
\[
 \boxed{(q_0,q_2,q_4,q_6,q_8)(F)
 =(0,-864,2016,0,0).}
 \tag{4.8}
\]
The second contraction moment is the following different linear
combination of the same five quadratic directions:
\[
 \boxed{\mu_2
 =126q_0+84q_2+36q_4+9q_6+q_8.}
 \tag{4.9}
\]
Thus its vanishing at \(F\) is the exact cancellation
\[
 84(-864)+36(2016)=0.
 \tag{4.10}
\]
By contrast,
\[
 I_2=\operatorname{tr}(A^2)
 =q_0+q_2+q_4+q_6+q_8,
 \tag{4.11}
\]
and (4.8) gives \(I_2(F)=1152\).

The invariant Hilbert function begins
\[
 \dim(R_4)_e=1,1,5,15,65,219
 \qquad(0\leq e\leq5).
 \tag{4.12}
\]
It follows that a minimal homogeneous generating set through degree three
has one generator in degree one, four new generators in degree two, and
ten new generators in degree three.  One may take \(\mu_1\) and any four
of the \(q_{2r}\) complementary to \(q_0\).  For completeness, the ten
primitive cubic multidegrees, indexed by the nontrivial summands in
(4.4), are
\[
\begin{split}
 &(112),(123),(134),(222),(224),\\
 &(233),(234),(244),(334),(444).
\end{split}
\tag{4.13}
\]
They are the unique trilinear contractions in those multidegrees.  With
the checker normalization by symmetrized operator traces, except for the
alternating ordered trace in multidegree \((234)\), their values at \(F\)
are
\[
\begin{array}{c|rrrrrrrrrr}
rst&112&123&134&222&224&233&234&244&334&444\\ \hline
c_{rst}(F)&-134784/5&497664/5&0&-10368&0&0&0&0&0&0 .
\end{array}
\tag{4.14}
\]
In particular, the cubic moment also misses nonzero cubic invariants, but
degree three is not minimal.

There is only one degree-one invariant, namely a nonzero scalar multiple
of \(\mu_1=\operatorname{tr}(A)\), and it vanishes at \(F\).  Equations
(4.8)--(4.11) therefore prove the sharp answer
\[
 \boxed{\min\{e:\exists I_e\in(R_4)_e,\ I_e(F)\ne0\}=2.}
 \tag{4.15}
\]
The coefficient determinant has
\[
 \det C=48,\qquad \det A=3981312,
 \tag{4.16}
\]
but it has degree five.  It is a separating invariant, not the first one.

### 4.2 The moment subalgebra

Let
\[
 {\cal A}=\mathbb Q[\mu_1,\mu_2,\ldots]\subset R_4.
 \tag{4.17}
\]
The low-degree comparison is already strict:
\[
\begin{array}{c|cc}
e&\dim{\cal A}_e&\dim(R_4)_e\\ \hline
1&1&1\\
2&2&5\\
3&3&15.
\end{array}
\tag{4.18}
\]
Here
\({\cal A}_2=\langle\mu_1^2,\mu_2\rangle\) and
\({\cal A}_3=\langle\mu_1^3,\mu_1\mu_2,\mu_3\rangle\).
Consequently the moments miss a three-dimensional quotient already in
degree two.  The class of \(I_2\) is one explicit missing direction, and
(4.10) describes exactly why it survives at \(F\).

There is also an exact global transcendence statement.  At the fixed
integer point stored by the checker, the \(22\) by \(25\) Jacobian of
\[
 \mu_1,\ldots,\mu_{22}
 \tag{4.19}
\]
has rank \(22\) modulo \(1000003\).  Hence these moments are algebraically
independent over \(\mathbb Q\).  Since
\(\dim R_4=25-3=22\), the moment algebra has full transcendence degree,
and the induced extension of fraction fields is finite algebraic.  In
terms of the presentation
\[
 \phi:\mathbb Q[t_1,t_2,\ldots]\longrightarrow R_4,\qquad
 t_m\longmapsto\mu_m,
 \tag{4.20}
\]
this proves
\[
 \ker\phi\cap\mathbb Q[t_1,\ldots,t_{22}]=0;
 \tag{4.21}
\]
relations must occur once further moments are adjoined, but their exact
ideal is not computed here.

The missing field direction is controlled by an equivariant involution.
Let \(S\) be the invariant apolar form on \(\operatorname{Sym}^4\),
written in the monomial basis as
\[
 S_{i,4-i}=1,-\frac14,\frac16,-\frac14,1
 \qquad(0\leq i\leq4),
 \tag{4.22}
\]
with all other entries zero, and define
\[
 \tau(A)=S^{-1}A^TS.
 \tag{4.23}
\]
This commutes with conjugation by \(\operatorname{SL}_2\).  In coefficient
coordinates it has the especially simple form
\[
 \boxed{\tau(C)_{ij}=(-1)^{i+j}C_{4-j,\,4-i}.}
 \tag{4.24}
\]
If \(C(x,y)=\sum c_{ij}x^iy^j\), then
\[
 \mu_m(C)=\sum_{I=0}^{4m}I!(4m-I)!
 [x^Iy^I]C(x,y)^m.
 \tag{4.25}
\]
Under (4.24), the diagonal coefficient indexed by \(I\) becomes the
coefficient indexed by \(4m-I\); its sign is
\((-1)^{2I}=1\), and the factorial weight is unchanged.  Therefore
\[
 \boxed{\mu_m\circ\tau=\mu_m\quad\hbox{for every }m\geq1.}
 \tag{4.26}
\]

On the five summands in (4.4), \(\tau\) acts respectively by
\[
 (+1,-1,+1,-1,+1).
 \tag{4.27}
\]
The primitive cubic \(c_{234}\) in (4.13) is consequently
\(\tau\)-odd, and the exact checker evaluates it nontrivially at a
rational test point.  Thus \(\tau\) acts nontrivially on \(R_4\), while
it fixes the entire moment algebra.  It follows that
\[
 \operatorname{Frac}{\cal A}
 \subseteq\operatorname{Frac}(R_4)^\tau
 \subsetneq\operatorname{Frac}R_4.
 \tag{4.28}
\]
In particular the generic moment map on the invariant quotient has degree
at least two: it cannot recover the apolar orientation.

At the special point \(F\), \(c_{234}(F)=0\).  This is structural rather
than an accidental extra cancellation: if
\[
 K=\operatorname{diag}(1,-1,1,-1,1)
 \tag{4.29}
\]
on \(\operatorname{Sym}^4\), exact evaluation gives
\[
 \tau(A_F)=K A_F K^{-1}.
 \tag{4.30}
\]
Thus \(\tau(F)\) is the
\(\operatorname{diag}(1,-1)\in\operatorname{PGL}_2\) translate of \(F\),
so every odd invariant vanishes there.  The degree-two separator
\(I_2\) is \(\tau\)-even.

Composing \(\tau\) with that translate gives a linear involution
\(\sigma\) fixing \(F\).  Its eigenspaces on \(V_4\) have dimensions
\[
 \dim V_4^+=15,\qquad \dim V_4^-=10.
 \tag{4.31}
\]
The three-dimensional orbit tangent splits as \(1+2\), so the
orbit-normal space \(V_4/\mathfrak {sl}_2F\) splits as
\[
 14+8.
 \tag{4.32}
\]
Every moment differential at \(F\) is \(\sigma\)-even.  The exact
all-order rank calculation in Section 3 becomes
\[
\begin{array}{c|cc}
&+&-\\ \hline
\text{rank of the moment differential}&12&0\\
\text{moment-fiber tangent modulo the orbit}&2&8.
\end{array}
\tag{4.33}
\]
Thus eight of the ten quotient-tangent directions in the all-moment fiber
are forced before any Hopf-specific calculation: pure moments cannot see
the apolar-odd normal directions.  Only the remaining two even directions
come from additional degeneracy of the moment map at this special point.

Despite the full transcendence degree, \(R_4\) is not integral over
\({\cal A}\).  Indeed, if the homogeneous positive-degree element \(I_2\)
were integral, graded integrality would give a homogeneous monic equation
for it whose remaining coefficients have positive degree in \({\cal A}\).
Evaluation at the common moment zero \(F\) would then give
\(I_2(F)^n=0\), contradicting (4.2).  The same argument applies to every
invariant nonzero at \(F\), including \(\det C\) and the three nonzero
cubic invariants in (4.14).  Thus the moment morphism is generically
algebraic but is not finite at the moment origin.

The conductor is now determined exactly:
\[
 ({\cal A}:R_4)=\{a\in{\cal A}:aR_4\subset{\cal A}\}
 =0.
 \tag{4.34}
\]
Indeed, if \(0\ne a\) belonged to the conductor, then for every
\(r\in R_4\) one would have \(ar\in{\cal A}\), hence
\(r=(ar)/a\in\operatorname{Frac}{\cal A}\), contradicting (4.28).
The exact relation ideal beyond (4.21) and the precise degree of the
generic moment map remain open.

## 5. First lifting obstructions

Let
\[
 X(\epsilon)=F+\epsilon H+\epsilon^2K+O(\epsilon^3).
\tag{5.1}
\]
For \(r\geq0\), the coefficient of \(\epsilon^2\) in the
\((r+1)\)-st moment, divided by \(r+1\), is
\[
 \ell_r(K)+q_r(H),\qquad
 q_0=0,\quad
 q_r(H)=\frac r2\mathcal E_2(H^2F^{r-1})\quad(r\geq1).
\tag{5.2}
\]
For each \(r\geq12\), let \(c_{rj}\) be the unique coefficients satisfying
\[
 \ell_r=\sum_{j=0}^{11}c_{rj}\ell_j.
\tag{5.3}
\]
After the first twelve equations have been used to eliminate \(K\), the
quadratic obstruction is
\[
 \boxed{
 \Omega_r(H)=q_r(H)-\sum_{j=0}^{11}c_{rj}q_j(H)=0.}
\tag{5.4}
\]
On the \(13\)-space (3.6), the seven quadrics
\[
 \boxed{\Omega_{12},\Omega_{13},\ldots,\Omega_{18}}
\tag{5.5}
\]
are linearly independent.  All later \(\Omega_r\) lie in their span; this
again follows by applying (3.3) to the bidegree-\((8,8)\) numerator
\(H^2\), not by declaring a cutoff.  Thus (5.5) are the complete first
(quadratic) lifting-obstruction equations.  Their exact coefficients in
the stored \(13\)-coordinate kernel basis are included in the generated
artifact; writing the large rational quadrics inline would obscure the
coordinate-free formula (5.4).

The orbit and radial directions satisfy all seven equations, as they must.

### 5.1 The reduced quadratic cone

Choose tangent coordinates so that the first four coordinates span the
three orbit directions and the radial direction.  Write
\(v_0,\ldots,v_8\) for the remaining quotient coordinates.  Exact
characteristic-zero primary decomposition of the seven quadrics gives a
scheme of dimension five and degree three.  Its radical is the linear
ideal
\[
\begin{aligned}
 v_0={}&0,\\
 35v_1-8v_4-70v_6+105v_7-105v_8={}&0,\\
 28v_2-43v_4+168v_5-336v_6+336v_7-336v_8={}&0,\\
 105v_3-251v_4+840v_5-1260v_6+1260v_7-1365v_8={}&0.
\end{aligned}
\tag{5.6}
\]
Thus the seven quadrics define a degree-three nonreduced thickening of one
five-plane.  In particular, the reduced quadratic tangent cone modulo
orbit and scaling is linear rather than a union of nonlinear branches.
The apolar involution gives this plane a sharper structure:
\[
 \boxed{\dim L^+=1,\qquad\dim L^-=4,}
 \tag{5.7}
\]
where \(L\) denotes the five-plane (5.6).  Restoring the even radial
direction, the reduced quadratic tangent cone modulo the orbit has split
\(2+4\).  Comparing with (4.33), the quadratic obstructions leave both
even moment-fiber directions and cut the eight odd directions down to
four.  This isolates the first nonlinear obstruction entirely inside the
apolar-odd excess.

### 5.2 Every reduced direction lifts through cubic order

Let \(H(h_0,\ldots,h_4)\) be the linear parametrization of (5.6) stored in
the artifact.  At the next order write
\[
 X(\epsilon)
 =F+\epsilon H+\epsilon^2K+\epsilon^3L+O(\epsilon^4).
\tag{5.8}
\]
After division by \(r+1\), the coefficient of \(\epsilon^3\) in moment
\(r+1\) is
\[
 \ell_r(L)
 +r\mathcal E_2(HKF^{r-1})
 +\frac{r(r-1)}6\mathcal E_2(H^3F^{r-2}).
\tag{5.9}
\]
Use the first twelve equations to choose a particular quadratic part of
\(K\).  Its remaining tangent part has thirteen coordinates.  Substitution
of a general quadratic polynomial in \(h_0,\ldots,h_4\) for those
coordinates turns (5.8) into a rational linear system with \(490\)
coefficient equations and \(195\) unknown coefficients.  Exact row
reduction has rank \(40\) and gives an explicit solution with \(40\)
nonzero coefficients.  Substituting that solution makes every cubic
obstruction identically zero.

The bidegree-\((12,12)\) version of the generating identity (3.3) supplies
the all-order tail reduction; orders \(12,\ldots,25\) are the exact
coefficient replay stored by the checker.  Hence
\[
\boxed{\text{every point of the reduced five-plane (5.6) lifts through
order }\epsilon^3.}
\tag{5.10}
\]
This includes every rank-drop boundary of the generic two-equation cubic
system; it is not only a generic-chart statement.

The \(1+4\) apolar split gives an immediate computational compression.
The tangent space in which the correction is chosen has split \(3+10\).
For a quadratic equivariant correction, the source monomials split as
\[
 \dim(\operatorname{Sym}^2L)^+=11,\qquad
 \dim(\operatorname{Sym}^2L)^-=4.
 \tag{5.11}
\]
Thus the \(195\) unrestricted correction coefficients reduce to
\[
 11\cdot3+4\cdot10=73.
 \tag{5.12}
\]
Only the eleven even cubic scalar monomials can occur in each of the
fourteen replay orders, reducing the existing \(490\) equations to
\(14\cdot11=154\).  The already-certified solution may be averaged under
\(\sigma\), so this smaller \(154\) by \(73\) system is guaranteed to be
consistent.

Consequently neither the quadratic nor cubic obstruction proves rigidity
modulo scaling.  The fourth-order coefficients contain \(K^2\), \(H^2K\),
\(HL\), and \(H^4\).  More precisely, for
\[
X(\epsilon)
=F+\epsilon H+\epsilon^2K+\epsilon^3L+\epsilon^4M+O(\epsilon^5),
\]
the coefficient of \(\epsilon^4\) in moment \(r+1\), divided by \(r+1\),
is
\[
\begin{aligned}
\ell_r(M)
{}&+r\mathcal E_2(HLF^{r-1})
  +\frac r2\mathcal E_2(K^2F^{r-1})\\
 &+\frac{r(r-1)}2\mathcal E_2(H^2KF^{r-2})
  +\frac{r(r-1)(r-2)}{24}\mathcal E_2(H^4F^{r-3}).
\end{aligned}
\tag{5.13}
\]
Equivariance reduces the cubic correction space to
\[
 11\cdot3+24\cdot10=273
 \tag{5.14}
\]
coefficients, because
\((\operatorname{Sym}^3L)^+\) and
\((\operatorname{Sym}^3L)^-\) have dimensions \(11\) and \(24\).
There are only
\[
 1+10+35=46
 \tag{5.15}
\]
even quartic monomials, according as they contain zero, two, or four odd
parameters.  The system therefore reduces canonically to at most
\[
 \boxed{644\text{ equations in }273\text{ correction coefficients}.}
 \tag{5.16}
\]

### 5.3 Fourth order: a section obstruction, but not an intrinsic one

There is a small all-order subtlety at fourth order.  Reducing the four
numerator degrees separately does not justify a cutoff, because their
indices and polynomial prefactors in (5.13) differ.  For arbitrary
balanced numerators \(A_d\) of degree \(d\), the relevant combined row is
\[
\begin{split}
 \Phi_r={}&\mathcal E_2(A_4F^r)
 +r\mathcal E_2(A_8F^{r-1})
 +\frac{r(r-1)}2\mathcal E_2(A_{12}F^{r-2})\\
&+\frac{r(r-1)(r-2)}{24}
  \mathcal E_2(A_{16}F^{r-3}).
\end{split}
\tag{5.17}
\]
Expanding the beta sum (3.3), now in degrees \(4,8,12,16\), and removing
the common nonzero hypergeometric factor puts every coordinate of
\(\Phi_r\) over one denominator with numerator degree at most \(41\).
Consequently
\[
 \boxed{\Phi_{12},\ldots,\Phi_{53}\text{ determine }\Phi_r
 \text{ for every }r\geq12.}
 \tag{5.18}
\]
The exact combined row rank is \(42\).  The checker uses the slightly
redundant range \(12,\ldots,57\) and replays the resulting rank through
order \(71\).  Thus the following conclusions concern the infinite moment
ideal, not a finite prefix.

First freeze the particular polynomial choices of \(K(h)\) and \(L(h)\)
made in Section 5.2, and seek a cubic polynomial tangent correction to
\(L\).  The resulting exact system has \(3220\) scalar coefficient
equations and \(455\) unknowns.  Its coefficient and augmented ranks are
\[
 90\quad\hbox{and}\quad91.
\tag{5.19}
\]
Thus this convenient global section of the lower-order lifting torsor
does not extend polynomially through fourth order.  On one generic
two-column chart its residual determinant is a product of two linear
forms and one quartic.  This is a **section obstruction** only: changing
the unused tangent part of the cubic lift changes it.

To test intrinsic liftability, take the nonzero reduced direction
\[
 (h_0,h_1,h_2,h_3,h_4)=(1,2,3,4,5)
 \tag{5.20}
\]
and restore all eleven free cubic-lift parameters left after the
third-order equations.  The fourth compatibility ideal has \(44\)
displayed generators.  Exact characteristic-zero Gröbner reduction
leaves one affine-linear equation and one quadric:
\[
 \boxed{\dim=9,\qquad\deg=2.}
 \tag{5.21}
\]
The quadric has rank one and discriminant square class \(41\).  Hence the
fiber is the union of two conjugate affine \(9\)-planes over
\(\mathbb Q(\sqrt {41})\); it has no \(\mathbb Q\)-point but is nonempty
over the algebraic closure.  In particular this nonradial reduced
direction has a fourth-order lift after a quadratic scalar extension.
Therefore fourth order still does not prove that the radial curve is the
complete reduced local quotient.

Equation (5.21) is only a four-jet statement.  It neither constructs a
formal arc nor proves that an additional reduced component exists; the
fifth and higher Kuranishi equations remain open.  Since
\(\mathcal E_2(ZF^m)\ne0\) for every \(m\geq1\), each of these mixed
contractions has nonzero constant term on any such infinitesimal lift and
is therefore a unit in its Artin base.  The lifted jets automatically
retain the \(Q=Z\) mixed defect.

The checker chooses one point on a \(\mathbb Q(\sqrt {41})\)-component
and reconstructs the complete jet
\[
 F+\epsilon H+\epsilon^2K+\epsilon^3L+\epsilon^4M.
 \tag{5.22}
\]
The supports of \(H,K,L,M\) in the \(25\)-monomial basis have sizes
\(16,13,11,9\).  Direct substitution verifies all four deformation
coefficients through the complete \(42\)-row tail.  Thus (5.21) is backed
by an explicit four-jet, not only a dimension computation.

### 5.4 The explicit nonradial four-jet is obstructed at fifth order

For a fifth correction \(N\), the coefficient of \(\epsilon^5\) in
moment \(r+1\), divided by \(r+1\), is
\[
\begin{aligned}
\ell_r(N)
{}&+r\mathcal E_2((HM+KL)F^{r-1})\\
&+\frac{r(r-1)}2\mathcal E_2((H^2L+HK^2)F^{r-2})\\
&+\frac{r(r-1)(r-2)}6\mathcal E_2(H^3KF^{r-3})\\
&+\frac{r(r-1)(r-2)(r-3)}{120}
  \mathcal E_2(H^5F^{r-4}).
\end{aligned}
\tag{5.23}
\]
The combined beta reduction now involves numerator degrees
\(4,8,12,16,20\).  After removing the common hypergeometric factor, its
common numerator degree is at most \(55\).  Hence the \(56\) rows
\[
 r=12,\ldots,67
 \tag{5.24}
\]
are all-order complete; their exact rank remains \(56\), with replay
through \(r=80\).

Restore the entire \(13\)-dimensional tangent freedom in \(M\).  On the
explicit four-jet (5.22), the fifth compatibility matrix and its
augmentation have ranks
\[
 \boxed{2\quad\hbox{and}\quad3.}
 \tag{5.25}
\]
There is a particularly short exact left obstruction.  In primitive
integer normalization its support on the residual rows
\((r=12,13,14)\) is
\[
 \boxed{(2727113757934325760,\,-407042494824,\,17047).}
 \tag{5.26}
\]
It annihilates all thirteen fourth-tangent columns and pairs nontrivially
with the constant fifth defect.  The pairing and its nonzero quadratic
norm are stored exactly.  Conjugation sends the certificate to the
conjugate explicit four-jet and leaves the norm nonzero.  Therefore
\[
\boxed{\text{both explicit conjugate nonradial four-jets fail at fifth
order}.}
\tag{5.27}
\]

This does not eliminate the two affine \(9\)-planes in (5.21): their
other points change \(K,L,M\), and the fifth obstruction may vanish on a
proper subvariety.  The next intrinsic calculation is the fifth
compatibility ideal over one whole \(\mathbb Q(\sqrt {41})\)-component,
followed by descent of its conjugate pair.

## 6. An exact positive-dimensional local quotient

For \(ab\ne0\), define
\[
\boxed{
 F_{a,b}=\frac12(aR+bZ)
 \left(2W(aR+bZ)^2-2abR^3-b^2R^2Z\right).}
\tag{6.1}
\]
On the unit sphere, with \(x=Z,\ t=T\), this becomes
\[
 p_{a,b}=\frac{a+bx}{2x}
 \left(a^2-t^2(a+bx)^2\right).
\tag{6.2}
\]
Put
\[
 J_{m,a}(X)=\int_0^X(a^2-v^2)^m\,dv.
\]
Then \(J'_{m,a}\) has a zero of order \(m\) at \(X=a\).  Coefficient
extraction at \(X=a+bx\) therefore proves, for every \(m\geq1\),
\[
\boxed{
\begin{aligned}
 \mathcal E_2(F_{a,b}^m)&=0,\\
 \mathcal E_2(ZF_{a,b}^m)
 &=(4m+2)!\,
   a^{2m+1}b^{m-1}\frac{m!}{(2m+1)!!}\ne0.
\end{aligned}}
\tag{6.3}
\]
This is the all-order generating identity, with no moment truncation.

The apparent two parameters split exactly into orbit and modulus.  If
\(t^{-2}=b/a\), then
\[
 \boxed{F_{a,b}=a^2b\,
 \bigl(\operatorname{diag}(t,t^{-1})\cdot F\bigr).}
\tag{6.4}
\]
Thus \(b/a\) is the diagonal \(G\)-orbit direction, while
\(\lambda=a^2b\) is radial.  Moreover
\[
 I_2(F_{a,b})=1152(a^2b)^2.
\tag{6.5}
\]
Consequently distinct \(\lambda\) near \(1\) are not
contraction-preserving equivalent.  Equations (6.3)--(6.5) construct a
one-dimensional family in the local quotient, every member retaining a
nonzero \(Q=Z\) mixed defect.

Hence the displayed counterexample is **not locally rigid**:
\[
\boxed{\dim_{\,[F]}({\cal Z}/G)\geq1.}
\tag{6.6}
\]
No claim is made here that the radial curve is the entire reduced local
quotient; deciding whether the additional obstructed tangent directions
support nonreduced or higher-order branches would require further local
obstruction calculations.

## Reproduction

Run

```bash
.venv/bin/python scripts/verify_two_pair_counterexample_local_moduli.py
.venv/bin/python scripts/verify_two_pair_counterexample_fourth_order.py
.venv/bin/python scripts/verify_two_pair_counterexample_fifth_order.py
```

The checker verifies the stabilizer, invariant, orbit tangent, exact
all-moment tangent rank and tail identities, the seven independent
quadratic obstruction equations, their degree-three five-plane radical,
the polynomial cubic lift of every reduced direction, the polynomial
family, and its pure and mixed identities.  It writes the exact matrices,
quadrics, and correction formulas to
`artifacts/generated-results/two_pair_counterexample_local_moduli.json`.
The second checker verifies the combined all-order fourth tail, the
obstructed polynomial section, and the nonradial fourth-order lift.  It
writes
`artifacts/generated-results/two_pair_counterexample_fourth_order.json`.
The third reconstructs an explicit \(\mathbb Q(\sqrt {41})\)-valued
four-jet and supplies its three-row all-order fifth obstruction in
`artifacts/generated-results/two_pair_counterexample_fifth_order.json`.
