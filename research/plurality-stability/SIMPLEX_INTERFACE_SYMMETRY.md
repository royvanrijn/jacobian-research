# Symmetry sectors of the Gaussian simplex interface

This note studies the full interface form from
[the positive-correlation continuation argument](POSITIVE_CORRELATION_CONTINUATION.md#the-full-interface-form-at-a-fixed-positive-correlation).
It gives the exact signed induction model for every cell count, identifies
the mass and rotation constraints, and proves positivity of the complete
three-cell sign sector for every positive correlation. An analytic
Hilbert--Schmidt reduction and rational interval certificate also give an
explicit full three-cell gap throughout $1/4\le\rho\le3/5$. The full
three-cell problem outside the known intervals, the higher-cell blocks,
and global optimality remain open. External review and formal verification
have not occurred.

## The signed facet representation

Identify the ambient space with $E=\mathbf1^\perp\subset\mathbb R^k$.
The simplex cells and facets are

\[
 T_i=\{x\in E:x_i\ge x_j\ (1\le j\le k)\},\qquad
 F_{ij}=\{x\in E:x_i=x_j\ge x_l\ (1\le l\le k)\}.
\]

Write $d\nu_{ij}=\phi_{k-1}\,d\mathcal H^{k-2}$ and extend facet
speeds antisymmetrically, $f_{ji}=-f_{ij}$. The Hilbert space, mass map,
and moment map are

\[
 \mathcal H=\bigoplus_{i<j}L^2(F_{ij},\nu_{ij}),\qquad
 (Cf)_i=\sum_{j\ne i}\int f_{ij}\,d\nu_{ij},\qquad
 (Df)_i=\sum_{j\ne i}\int x f_{ij}(x)\,d\nu_{ij}(x).
 \tag{S1}
\]

Thus $\mathcal H_0=\ker C$ is the balanced space. The permutation action is

\[
 (U_gf)_{ij}(x)=f_{g^{-1}i,g^{-1}j}(g^{-1}x).
 \tag{S2}
\]

It is unitary, and $C,D$ are equivariant. Permutations preserve the
Mehler kernel and carry $a_{\rho,ij}(x)$ to
$a_{\rho,gi,gj}(gx)$. Consequently the full form $\mathfrak q_\rho$
and its compression to $\mathcal H_0$ commute with $S_k$ for every
$0<\rho<1$.

**Signed induction formula.** For $k\ge3$, let
$K=L^2(F_{12},\nu_{12})$ with the geometric $S_{k-2}$ action on the
remaining coordinates. Over $\mathbb C$,

\[
 \mathcal H_{\mathbb C}\simeq
 \operatorname{Ind}_{S_2\times S_{k-2}}^{S_k}
       (\operatorname{sgn}_2\boxtimes K_{\mathbb C}).
 \tag{S3}
\]

Indeed the stabilizer of the unoriented facet $\{1,2\}$ is
$S_2\times S_{k-2}$. The transposition $(12)$ fixes its supporting
plane pointwise but reverses the speed orientation. It acts by minus
the identity. Representatives of the stabilizer cosets carry the
reference facet to all facets, giving (S3), with the invariant surface
measures supplying the Hilbert-space identification. Treating $(12)$ as
acting trivially would give the wrong decomposition.

Decompose $K_{\mathbb C}=\bigoplus_{\mu\vdash k-2}V_\mu\otimes M_\mu$.
Then

\[
 \mathcal H_{\mathbb C}\simeq
 \bigoplus_{\lambda\vdash k}V_\lambda\otimes
 \left(\bigoplus_{\substack{\mu\vdash k-2\\
                   \lambda/\mu\text{ a vertical two-strip}}}M_\mu\right).
 \tag{S4}
\]

The algebraic inputs are the characteristic map for induced characters
and the vertical Pieri rule: Macdonald, *Symmetric Functions and Hall
Polynomials*, second edition, Chapter I, (7.3), (7.5)--(7.6), printed
pp. 113--114, and (5.17), printed pp. 73--74. The theorem statements and
the alternant proof of (5.17) were inspected in the
[primary text](https://math.berkeley.edu/~corteel/MATH249/macdonald.pdf).
Applying them to $e_2s_\mu$ gives (S4).

Every $M_\mu$ is infinite dimensional. Choose a small open ball in the
relative interior of $F_{12}$ where the other coordinates are distinct,
with disjoint translates under $S_{k-2}$. Functions on the union of these
balls carry the regular representation tensored with the infinite
dimensional space of functions on one ball. Separability gives countable
multiplicity. Every partition $\lambda$ with at least two rows admits a
vertical two-strip removal: remove the last box of each of its bottom
two rows. The result is still a partition. The one-row partition admits
no such removal. Therefore **every nontrivial irreducible representation
occurs with infinite multiplicity, and the trivial representation is
absent**. For $k=2$, the facet is a point and $\mathcal H_0=0$ instead.

Schur decomposition thus gives finitely many representation types, each
with an operator on an infinite dimensional multiplicity space. It does
not reduce the spectral problem to finitely many scalar eigenvalues.

## The constrained copies and first-moment information

**Mass and rotation copies.** Put $v_k=\nu_{ij}(F_{ij})>0$, independent
of the pair. For $p\in E$,

\[
 (C^*p)_{ij}=p_i-p_j,\qquad CC^*p=kv_kp.
 \tag{S5}
\]

Hence $C$ removes exactly one standard copy $V_{(k-1,1)}$:
$\mathcal H_0=\mathcal H\ominus C^*E$. The geometric rotation space is

\[
 \mathcal R_k\simeq\Lambda^2E\simeq V_{(k-2,1,1)}.
 \tag{S6}
\]

For the last identity, the constant signed-edge module is
$\Lambda^2\mathbb C^k=\operatorname{Ind}(\operatorname{sgn}_2\boxtimes
\mathbf1)$, whose vertical Pieri decomposition is
$V_{(k-1,1)}\oplus V_{(k-2,1,1)}$. Also
$\Lambda^2(\mathbf1\oplus E)=E\oplus\Lambda^2E$.
The map from a skew operator $\Omega$ to
$f^\Omega_{ij}(x)=-(\Omega x)\cdot n_{ij}$ is injective:
if it vanishes on every facet interior, each root normal is an eigenvector
of $\Omega$; skew symmetry makes its eigenvalue zero, and the roots span
$E$. This proves the geometric identification. When $k=3$ the hook in
(S6) is the sign representation. Removing rotations removes one copy,
not its entire isotypic space.

For $k\ge4$ the moment target has the decomposition

\[
 E\otimes E=\mathbf1\oplus E\oplus V_{(k-2,2)}
                         \oplus V_{(k-2,1,1)}.
 \tag{S7}
\]

One derivation uses the horizontal Pieri decomposition of the permutation
module on unordered pairs, together with
$\operatorname{Sym}^2\mathbb C^k=\mathbb C^k\oplus
\mathbb C^{\binom k2}$. For $k=3$, omit $V_{(k-2,2)}$.
Since the domain has no trivial component, $D$ vanishes on all sectors
other than the three nontrivial types in (S7). Even in those types it
sees only finite dimensional subspaces of their multiplicity spaces.

**Exact rotation gauge.** The
[first-moment interface theorem](SMALL_POSITIVE_CORRELATION_INTERVAL.md#extension-to-every-fixed-cell-count)
gives $Q_0(f)=\ell\|f\|^2-\|Df\|_F^2\ge0$ on $\mathcal H_0$
with kernel $\mathcal R_k$. Polarization at $r\in\mathcal R_k$ gives

\[
 \langle Dr,Df\rangle_F=\ell\langle r,f\rangle_{\mathcal H}.
 \tag{S8}
\]

As $D:\mathcal R_k\to\mathcal V$ is an isomorphism onto the moment
rotation space, the moment gauge is exactly

\[
 \{f\in\mathcal H_0:Df\perp\mathcal V\}
             =\mathcal H_0\cap\mathcal R_k^\perp.
 \tag{S9}
\]

Thus the nonrotational copies of the hook are invisible to $D$ as well.
On every first-moment-invisible space $Q_0=\ell I$. This statement
concerns the first-moment form. Higher Hermite degrees still contribute
to $\mathfrak q_\rho$ at positive correlation.

For four cells, let $K_+$ and $K_-$ be the even and odd functions on
$F_{12}$ under $(34)$. Formula (S4) reads

| $S_4$ type | Multiplicity space | Constraint or known zero mode |
|---|---|---|
| $(3,1)$ | $K_+$ | Remove one mass copy |
| $(2,2)$ | $K_-$ | No geometric zero mode |
| $(2,1,1)$ | $K_+\oplus K_-$ | Remove one rotation copy |
| $(1,1,1,1)$ | $K_-$ | First-moment invisible |

Each multiplicity space is infinite dimensional. The full operator may
mix the two branches inside the $(2,1,1)$ multiplicity space. Apart from
the removed finite dimensional copies, all four types remain possible
locations of an additional zero mode.

## Interaction sparsity and the sign-sector Hermite floor

Orient an edge $e=(i,j)$ and put $d_e=e_i-e_j$ in label space. Direct
expansion of the measures $\lambda_i^f$ gives

\[
 B_\rho(f)=\frac1\rho\sum_{e,e'}(d_e\cdot d_{e'})
 \iint_{F_e\times F_{e'}}M_\rho(x,y)f_e(x)f_{e'}(y)
                                  \,d\nu_e(x)d\nu_{e'}(y).
 \tag{S10}
\]

The incidence coefficient is $2$ on the diagonal, $\pm1$ for edges
sharing a vertex, and zero for disjoint edges. This is an exact sparsity
rule for the interaction; it does not make an individual symmetry block
finite dimensional.

Let $H_\alpha$ be orthonormal Gaussian Hermite polynomials. Mehler
expansion yields

\[
 B_\rho(f)=\sum_{m\ge0}\rho^{m-1}E_m(f),\qquad
 E_m(f)=\sum_i\sum_{|\alpha|=m}
                \left|\int H_\alpha\,d\lambda_i^f\right|^2\ge0.
 \tag{S11}
\]

To justify the expansion for surface measures, the absolute Hermite
series at $(x,y)$ is at most
$\sqrt{M_\rho(x,x)M_\rho(y,y)}$. Its integral against the total
variation measures is finite by Cauchy--Schwarz on the facets:
$M_\rho(x,x)$ grows as
$\exp(\rho|x|^2/(1+\rho))$, whose exponent is below $|x|^2/2$.
This permits termwise integration. Balance gives $E_0=0$.

**Sign-sector floor.** For $k\ge3$, on the entire sign isotypic space,

\[
 E_m=0\quad\text{if }m<d_k:=\frac{(k-1)(k-2)}2.
 \tag{S12}
\]

Indeed the coefficient map of degree $m$ is equivariant into
$E\otimes\operatorname{Sym}^mE$. A sign vector in this target is a
centered tuple of polynomials $P_i$ with
$P(gx)=\operatorname{sgn}(g)gP(x)$. Its first component must be
alternating in $x_2,\ldots,x_k$. On $E$ these are independent
coordinates, so it is divisible by
$\prod_{2\le i<j\le k}(x_i-x_j)$, of degree $d_k$. This proves (S12).

The degree is attained. With a consistent Vandermonde convention take
$P_i=(-1)^{i-1}\Delta(x_1,\ldots,\widehat{x_i},\ldots,x_k)$.
The alternating-minor identity gives $\sum_iP_i=0$; it is the determinant
with two constant rows. The tuple has the required transformation law.
Its lower Hermite components vanish by the same degree obstruction.
Set $f_{ij}=P_i-P_j$ on $F_{ij}$. This is a sign speed and is balanced,
because $C$ maps into the standard representation. It is nonzero: on
$F_{12}$, $P_2=-P_1$ and $P_1$ is not identically zero. Pairing its
degree-$d_k$ coefficient with $P$ gives
$\sum_{i<j}\int(P_i-P_j)^2d\nu_{ij}>0$.

In particular, for $0<\rho\le b<1$ and sign speeds,

\[
 0\le B_\rho(f)\le(\rho/b)^{d_k-1}B_b(f).
 \tag{S13}
\]

For four cells the interaction starts at normalized order $\rho^2$;
for five cells at order $\rho^5$. No arbitrary-correlation sign conclusion
follows from this suppression alone.

## Three cells: a full positive-correlation sign-sector theorem

Now $k=3$ and the facets are rays. Orient them cyclically as $(12)$,
$(23)$, $(31)$, and let $r\ge0$ be distance from the junction. Their
unit tangents have pairwise inner product $-1/2$. Write

\[
 d\nu(r)=(2\pi)^{-1}e^{-r^2/2}\,dr,\quad
 \sigma=\sqrt{1-\rho^2},\quad c=\rho/\sigma,\quad\phi=\phi_1,
\]

\[
 M_t(r,s)=\sigma^{-2}\exp\left(
  -\frac{\rho^2(r^2+s^2)}{2\sigma^2}
  +\frac{\rho t rs}{\sigma^2}\right),\qquad
 K_-=2(M_1-M_{-1/2}),\quad K_+=2M_1+M_{-1/2}.
 \tag{S14}
\]

The sign sector consists of $f=(h,h,h)/\sqrt3$ and is automatically
balanced. The standard sector consists of sum-zero edge vectors; balance
requires zero $\nu$-mean in each of its two scalar components. The
incidence matrix in (S10) is diagonal $2M_1$, off-diagonal
$-M_{-1/2}$. Thus its sign and standard scalar kernels are exactly
$K_-$ and $K_+$.

The common local coefficient is

\[
 a_\rho(r)=\frac\rho\sigma\left[
  2\phi(0)\Phi(cr)+\phi(\sqrt3cr/2)\Phi(-cr/2)\right].
 \tag{S15}
\]

To derive it, use inward-normal and ray-tangent coordinates $(t,z)$,
where $T_1=\{t\ge0,z\ge-t/\sqrt3\}$. The reflected derivative is
$2\rho\sigma^{-1}\int_0^\infty t\phi(t)\Phi(cr+t/\sqrt3)dt$.
Integration by parts and completion of the square give (S15).
In particular $a_\rho/\rho\ge\ell/\sigma$ with
$\ell=3\phi(0)/2$.

Completion of the square also gives, for $t=1,-1/2$,

\[
 \int M_t(r,s)\,d\nu(s)
  =\frac{\phi(0)}\sigma e^{-(1-t^2)c^2r^2/2}\Phi(tcr),
\]

\[
 \int M_t(r,s)s\,d\nu(s)
  =\frac{\phi(0)}\sigma e^{-(1-t^2)c^2r^2/2}
       [\rho tr\Phi(tcr)+\sigma\phi(tcr)].
 \tag{S16}
\]

The identity
$e^{-3u^2/8}\phi(u/2)=\phi(u)$ then proves

\[
 \int K_-(r,s)s\,d\nu(s)=a_\rho(r)r,\qquad
 \int K_+(r,s)\,d\nu(s)=a_\rho(r)/\rho.
 \tag{S17}
\]

**All-correlation sign gap.** For every $0<\rho<1$, the sign-sector
form is nonnegative, its kernel is exactly the one-dimensional geometric
rotation space, and it has a positive gap modulo that space.

For all $h\in L^2(\nu)$, the first identity in (S17) gives

\[
 \rho\,\mathfrak q_\rho((h,h,h)/\sqrt3)
 =\frac12\iint K_-(r,s)rs
       \left(\frac{h(r)}r-\frac{h(s)}s\right)^2d\nu(r)d\nu(s).
 \tag{S18}
\]

The terms are integrable: each diagonal term is
$\int a_\rho h^2d\nu<\infty$ by (S17), and Cauchy--Schwarz bounds
the cross term. The value at $r=0$ is irrelevant. The kernel $K_-$ is
strictly positive for $r,s>0$, so equality forces $h(r)/r$ to be constant
almost everywhere. These are precisely rotational speeds. The positive
gap follows from the positive multiplication coefficient and compact
interaction, as in (B8) of the continuation note. No maximizing assumption
is used in this proof.

There is also an explicit estimate beyond the perturbative argument.
Let $T=a_\rho^{-1/2}K_-a_\rho^{-1/2}$ on $L^2(\nu)$. It is compact,
nonnegative and has simple top eigenvalue $1$, by (S18); its eigenvector
is $\sqrt{a_\rho(r)}r$. Nonnegativity also follows by expanding the
kernel: every coefficient of
$2(e^z-e^{-z/2})$ is nonnegative for $z\ge0$. Define

\[
 J_\rho(u)=\frac{\frac14+\frac1{2\pi}
                   \arcsin(\rho u/(1+\rho^2))}
                 {\sqrt{(1+\rho^2)^2-\rho^2u^2}},\qquad
 H(\rho)=\frac{16}{9\rho^2}
       [J_\rho(2)-2J_\rho(1/2)+J_\rho(-1)].
 \tag{S19}
\]

Integrating a two-variable Gaussian on a quadrant and using
$a_\rho\ge\rho\ell/\sigma$ gives $\|T\|_{HS}^2\le H(\rho)$.
More explicitly,
$\iint M_tM_u\,d\nu d\nu=J_\rho(t+u)/(2\pi\sigma^2)$.
If $H(\rho)<2$, the second eigenvalue is at most
$\sqrt{H(\rho)-1}$, since the top eigenvalue contributes $1$ to the
Hilbert--Schmidt square. Changing back from $\sqrt{a_\rho}h$ gives

\[
 \mathfrak q_\rho((h,h,h)/\sqrt3)
 \ge \frac\ell\sigma(1-\sqrt{H(\rho)-1})
                    \inf_{b\in\mathbb R}\|h-br\|_{L^2(\nu)}^2.
 \tag{S20}
\]

At $\rho=1/2$ the entirely elementary bound $H(1/2)<7/4$ holds.
Indeed $3<\pi<22/7$, $\arcsin(4/5)<15/16$ (use
$\sin t\ge t-t^3/6$), and $\arcsin x\ge x$ give

\[
 H(1/2)\le\frac{64}9\left[
  \frac{13}{24}-\frac{62}{55\sqrt6}+\frac{41}{55\sqrt{21}}\right]
 <\frac{64}9\left[\frac{13}{24}-\frac{1240}{2695}
                                  +\frac{451}{2750}\right]<\frac74.
 \tag{S21}
\]

Here $1/\sqrt6>20/49$ and $1/\sqrt{21}<11/50$; the final
comparison is rational arithmetic. Consequently

\[
 \boxed{\displaystyle
 \mathfrak q_{1/2}((h,h,h)/\sqrt3)
 \ge\ell(2/\sqrt3-1)\inf_b\|h-br\|_{L^2(\nu)}^2.}
 \tag{S22}
\]

This is a sector gap, not a full-interface or global-optimality theorem.

## The standard sector and an explicit full gap

The full three-cell interface-kernel assertion is now equivalent to the
following strict scalar inequality for every $0<\rho<1$:

\[
 \int a_\rho(r)h(r)^2d\nu(r)
   >\iint K_+(r,s)h(r)h(s)d\nu(r)d\nu(s),\qquad
 0\ne h\in L^2(\nu),\quad\int h\,d\nu=0.
 \tag{S23}
\]

For fixed $\rho$, strict positivity implies a positive gap by compactness
and the positive lower bound for $a_\rho$. The sign sector has already
been settled by (S18). The second identity in (S17) describes the
unbalanced translation mode, with generalized eigenvalue $1/\rho$.
The mass constraint uses $\nu$-mean, whereas its natural spectral
orthogonality uses $a_\rho\nu$-mean. Removing the top generalized
eigenfunction without accounting for this difference does not prove
(S23). Also $K_+-K_-=3M_{-1/2}$ is not asserted to be a positive or
negative quadratic form on this mean-zero space.

There is nevertheless a direct bound on a substantial explicit interval.
Put $v=\nu([0,\infty))=\phi(0)/2$, and define

\[
 P(x)=\frac14+\frac{\arcsin x}{2\pi},\quad
 R_\rho=4J_\rho(2)+4J_\rho(1/2)+J_\rho(-1),\quad
 L_\rho=\frac{2P(\rho)}{\sqrt{1-\rho^2}}
                 +\frac{P(-\rho/2)}{\sqrt{1-\rho^2/4}},
\]

\[
 G(\rho)=\frac4{9\rho^2}R_\rho
              -\frac{16(1-\rho^2)}{9\rho^2}L_\rho^2.
 \tag{S24}
\]

Let $c_\rho=v^{-2}\iint K_+\,d\nu d\nu$ and
$\overline K_+=K_+-c_\rho$. On mean-zero $h$ subtracting this constant
does not change the quadratic interaction. The same quadrant Gaussian
integral as in (S19) gives

\[
 \iint K_+\,d\nu d\nu=\frac{L_\rho}{2\pi},\qquad
 \iint K_+^2\,d\nu d\nu=\frac{R_\rho}{2\pi(1-\rho^2)},
 \qquad
 \frac{\|\overline K_+\|_{HS}^2}{(\rho\ell/\sigma)^2}=G(\rho).
 \tag{S25}
\]

Here the squared norm after subtracting the optimal constant is
$\|K_+\|_{HS}^2-v^{-2}(\iint K_+)^2$. These formulas can also be
derived in polar coordinates from
$\int_0^{\pi/2}(A-B\sin2\theta)^{-1}d\theta$.
Thus, whenever $G(\rho)<1$, the standard sector satisfies

\[
 \mathfrak q_\rho(h)\ge\frac\ell\sigma(1-\sqrt{G(\rho)})
                                      \|h\|_{L^2(\nu)}^2.
 \tag{S26}
\]

This uses the ordinary $L^2(\nu)$ operator norm of the centered kernel,
bounded by its Hilbert--Schmidt norm; it does not replace the mass
constraint by weighted spectral orthogonality.

**Explicit three-cell full gap.** For every $1/4\le\rho\le3/5$ and
every balanced $f\in\mathcal H_0$,

\[
 \boxed{\displaystyle
 \mathfrak q_\rho(f)\ge
  \ell(1-3/\sqrt{10})\operatorname{dist}(f,\mathcal R_3)^2.}
 \tag{S27}
\]

In particular rotations are the only zero modes throughout this interval.
At $\rho=1/2$ the stronger bound is

\[
 \mathfrak q_{1/2}(f)\ge
  \ell(2/\sqrt3-1)\operatorname{dist}(f,\mathcal R_3)^2.
 \tag{S28}
\]

The exact certificate verifies $H(\rho)<19/10$ and $G(\rho)<4/5$
throughout $[1/4,3/5]$. Equations (S20) and (S26), followed by the
orthogonal sign/standard decomposition, give (S27), since
$\sqrt{4/5}<3/\sqrt{10}$ and $\sigma\le1$. At $1/2$, the certificate
also gives $G(1/2)<9/16$; combine (S22) with (S26) to get (S28).

The certificate is a finite reduction covering the whole stated interval,
not evaluations at sample correlations. It uses the 350 closed intervals
$[j/1000,(j+1)/1000]$, $250\le j<600$, and rational interval arithmetic
in the displayed formulas for $H$ and $G$. Square roots are enclosed by
integer squares. For $0\le x<1$ the positive series

\[
 \arcsin x=\sum_{n\ge0}\frac{\binom{2n}n}{4^n(2n+1)}x^{2n+1}
\]

has decreasing coefficients; after 64 terms its tail is bounded by the
next coefficient times $x^{129}/(1-x^2)$. Odd symmetry handles negative
arguments. Alternating arctangent sums and
$\pi=16\arctan(1/5)-4\arctan(1/239)$ enclose $\pi$. The identity
follows from $\tan(4\arctan(1/5))=120/119$ and the subtraction formula.
Every subsequent rounding is outward to rational endpoints. The checker
and receipt below retain the exact grid, bounds, algorithm and hashes;
replaying this calculation is part of verifying (S27)--(S28).

For $k\ge4$, (S4) supplies the precise symmetry blocks but no sign for
their remaining infinite dimensional parts. In particular (S18) does
not extend merely by replacing the three rays by higher dimensional
facets: a full sign-type function may change sign inside a reference
facet.

## Monotonicity and the other global obligations

**Raw-deficit obstruction.** For every $k\ge2$, the unnormalized
simplex deficit $\Delta_\rho(A)=S_\rho(T)-S_\rho(A)$ need not be
monotone on $(0,1)$. Take a balanced radial partition in $\mathbb R^{k-1}$
using successive radius quantiles. All its first moments vanish, whereas
$W_1(T)=F_k>0$. The Hermite expansion gives
$\Delta_\rho(A)=F_k\rho+O(\rho^2)>0$ for small positive $\rho$.
At $\rho=0$ both stability values are $1/k$. Strong continuity of the
Ornstein--Uhlenbeck semigroup at $\rho=1$ gives both values tending to
$1$. Thus the deficit tends to zero at both endpoints and is positive
somewhere between them, ruling out monotonicity in either direction.
For $k=2$ one may use the centered interval of Gaussian measure $1/2$
against the halfline; then $F_2=1/\pi$.

A useful monotonicity principle therefore needs a specified normalization
or a different quantity. The perimeter-scale candidate
$\Delta_\rho/\arccos\rho$ is not proved monotone here.

The remaining targets are distinct:

* Prove (S23), and obtain explicit bounds in all the blocks of (S4).
* Classify stationary stable balanced partitions. Simplex spectral
  stability does not classify critical points away from the simplex.
* Exclude a nonsimplex maximizer tied with the simplex. The continuation
  dichotomy retains this alternative even if the entire simplex Hessian
  is strictly positive modulo rotations.
* Prove exact optimality near $\rho=1$ using a quantitative perimeter
  bound together with estimates relative to the deformation of a moving
  competitor. A perimeter limit for each fixed partition is insufficient.
* Seek stronger moment or perimeter inequalities that control that
  relative error, and formulate a normalized monotonicity statement
  before testing it.

Run the exact checker and its optional exploratory quadrature with

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python research/plurality-stability/check_simplex_interface_symmetry.py --quadrature --output research/plurality-stability/simplex_interface_symmetry_verification.json
```

The exact part requires only Python's standard library. It checks all 30
permutations in the $S_3/S_4$ character examples, (S21), and the interval
inequalities used in (S27)--(S28), including rejection controls. The
analytic operator reduction and the all-$k$ representation argument
still require written proof review. The optional quadrature requires
NumPy 1.26.4 and evaluates the open standard-sector inequality at four
chosen correlations with two grids on $[0,10]$. It has no certified
tail or discretization errors and proves no inequality. Its agreement
is separate from the rational interval certificate.
