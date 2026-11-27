# Symmetry sectors of the Gaussian simplex interface

This note studies the full interface form from
[the positive-correlation continuation argument](POSITIVE_CORRELATION_CONTINUATION.md#the-full-interface-form-at-a-fixed-positive-correlation).
It gives the exact signed induction model for every cell count, identifies
the mass and rotation constraints, and proves positivity of the complete
three-cell sign sector for every positive correlation. Reflected
one-dimensional Ornstein--Uhlenbeck comparisons and a single elementary
bound at correlation $1/2$ now give an explicit full three-cell gap for
**every $0<\rho<1$**, uniformly in the ambient dimension $d\ge2$.
The complete proof is in the section on reflected kernels below. For every
cell count, the ambient-dimension extension reduces exactly to the balanced
problem in the active space, with an explicit loss in the gap. The higher-cell
sign and standard kernels also have exact positive determinant and
permanent formulas on an ordered chamber. In the rotation hook, positivity
in chamber-gap coordinates and a left/right ground-state argument prove
the rotational-only kernel and a positive gap at every correlation for
every $k$. Strict supersolutions also settle the $(k-2,2)$ and sign
sectors for all $k\ge4$, with the explicit sign bound (S62) for $k\ge5$.
The exterior Gaussian identity (S72) and a linear positive
supersolution also settle the whole third exterior hook
$(k-3,1,1,1)$ for every $k\ge4$, by (S73).
The four-cell problem reduces
to one scalar standard-sector inequality. The exact polynomial
certificate (S74)--(S76) also settles the five-cell $(2,2,1)$ block,
so five cells reduce to the standard inequality as well.
These results do not supply
the missing full active-space gap
for $k\ge4$ at arbitrary correlation. The remaining higher-cell blocks and
global optimality remain open. External review and formal verification
have not occurred.
The common compressed-immanant formula (S67)--(S69) gives entrywise
nonnegative matrix kernels in a fixed Young basis for every cell count
and every symmetry block. It leaves the remaining spectral bounds explicit.

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
| $(2,2)$ | $K_-$ | Strict gap by (S60) |
| $(2,1,1)$ | $K_+\oplus K_-$ | Rotation copy; gap off it by (S56) |
| $(1,1,1,1)$ | $K_-$ | First-moment invisible; strict gap by (S61) |

Each multiplicity space is infinite dimensional. The full operator may
mix the two branches inside the $(2,1,1)$ multiplicity space. The complete
hook argument (S52)--(S56) below handles that mixing and excludes every
additional hook zero mode at every positive correlation. The $(2,2)$
and sign types are strictly positive by (S60)--(S61). Only the standard
type remains to be estimated.

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
for five cells at order $\rho^5$. The all-correlation sign gap is proved
in (S62)--(S66), using additional total-positivity and Gaussian-matrix
arguments; suppression alone would not prove it.

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

The full three-cell interface-kernel assertion in the active plane is
equivalent to the following strict scalar inequality, proved below for
every $0<\rho<1$:

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
replaying this calculation is part of verifying this particular proof
of (S27)--(S28). The all-correlation proof below needs only the bound
$G(1/2)<9/16$, for which an elementary certificate is supplied separately.

## Reflected kernels close the three-cell spectral problem

**Boundary interaction monotonicity.** For every fixed cell count and
balanced speed, (S11) gives

\[
 B_\rho(f)\le B_b(f)\qquad(0<\rho\le b<1).
 \tag{S29}
\]

This is monotonicity of the normalized boundary interaction, not of the
noise-stability deficit. All coefficients in (S11) are nonnegative and
the degree-zero coefficient vanishes by balance. The convergent series
therefore proves (S29) on the entire balanced Hilbert space.

For three cells introduce the two reflected kernels

\[
 N=M_1+M_{-1},\qquad D=M_1-M_{-1}.
\]

Here $M_{-1}$ is the same explicit formula (S14), evaluated at $t=-1$;
it represents reflection to the opposite halfline, not a simplex facet.
The even and odd extensions of $h$ from $[0,\infty)$ to $\mathbb R$
identify their operators with $\phi(0)/\sigma$ times the ordinary
one-dimensional OU operator. Consequently

\[
 \int N(r,s)d\nu(s)=\frac{\phi(0)}\sigma,\qquad
 \int D(r,s)s\,d\nu(s)=\frac{\rho\phi(0)}\sigma r.
 \tag{S30}
\]

For clarity, the ordinary OU spectrum used here follows from the
generating function $e^{tx-t^2/2}=\sum_mH_m(x)t^m/\sqrt{m!}$:
$T_\rho e^{tx-t^2/2}=e^{\rho tx-\rho^2t^2/2}$, so
$T_\rho H_m=\rho^mH_m$. Gaussian integration of two generating
functions gives orthonormality. If an $L^2(\gamma_1)$ function is
orthogonal to every $H_m$, its Gaussian Laplace transform is entire
with all Taylor coefficients zero; evaluation on imaginary arguments
and uniqueness of Fourier transforms give that the function is zero.
This proves completeness as well. Gaussian exponential integrability
and Cauchy--Schwarz justify these operations.

For even extensions of zero mean, the first allowed Hermite degree is
$2$. For odd extensions orthogonal to $r$, it is $3$. The OU eigenvalues
$\rho^m$ therefore give

\[
 \frac12\iint N(r,s)(h(r)-h(s))^2d\nu d\nu
 \ge\frac{\phi(0)}\sigma(1-\rho^2)
                  \inf_b\|h-b\|^2,
\]

\[
 \frac12\iint D(r,s)rs
       (h(r)/r-h(s)/s)^2d\nu d\nu
 \ge\frac{\phi(0)}\sigma(\rho-\rho^3)
                  \inf_b\|h-br\|^2.
 \tag{S31}
\]

The norms are in $L^2(\nu)$. These are Hilbert-space statements for all
$h\in L^2(\nu)$, without trace or differentiability assumptions.

The decisive pointwise comparisons are

\[
 K_+\ge\tfrac32N,\qquad K_-\ge\tfrac32D
                       \qquad(r,s\ge0).
 \tag{S32}
\]

For $z=\rho rs/\sigma^2$ and $w=e^{z/2}\ge1$, multiply by the
common positive Gaussian factor and by $2e^z$. The two differences
become respectively

\[
 w^4+2w-3=(w-1)(w^3+w^2+w+3),\qquad
 w^4-4w+3=(w-1)^2(w^2+2w+3).
\]

Thus (S32) is pointwise positivity, not a claimed Loewner comparison of
the kernels themselves. It can be integrated against the nonnegative
squared differences in (S18) and (S31).

For the sign sector, (S18), (S31) and (S32) immediately give the explicit
all-correlation bound

\[
 \mathfrak q_\rho((h,h,h)/\sqrt3)
       \ge\ell\sigma\inf_b\|h-br\|^2.
 \tag{S33}
\]

For the standard sector put $R(r)=\int K_+(r,s)d\nu(s)=a_\rho(r)/\rho$.
It is positive and at most $2\phi(0)/\sigma$. The exact rewrite is

\[
 \rho\mathfrak q_\rho(h)=
 \frac12\iint K_+(r,s)(h(r)-h(s))^2d\nu d\nu
                  -(1-\rho)\int R h^2d\nu.
\]

For $\int h\,d\nu=0$, use (S31)--(S32) and the upper bound on $R$:

\[
 \mathfrak q_\rho(h)\ge
   \frac{\phi(0)(1-\rho)(3\rho-1)}{2\rho\sigma}\|h\|^2
   =\frac{\ell(1-\rho)(3\rho-1)}{3\rho\sigma}\|h\|^2.
 \tag{S34}
\]

This is positive for every $\rho>1/3$, including arbitrarily near one.
For $0<\rho\le1/2$, apply (S29) in the balanced standard sector and
use the centered-kernel bound at $b=1/2$:

\[
 B_\rho(h)\le B_{1/2}(h)
 \le\frac{\ell}{\sqrt{3}/2}\sqrt{G(1/2)}\|h\|^2
 <\frac{\sqrt3}2\ell\|h\|^2.
\]

The lower bound $A_\rho\ge\ell\sigma^{-1}I$ proves

\[
 \mathfrak q_\rho(h)\ge
  \ell(\sigma^{-1}-\sqrt3/2)\|h\|^2\qquad(0<\rho\le1/2).
 \tag{S35}
\]

The two ranges overlap; there is no uncovered correlation interval.

**An elementary certificate at $1/2$.** This proof does not require the
350-interval calculation. Directly from (S24), write
$G(1/2)=16\mathscr R/9-16\mathscr L^2/3$, where

\[
 \mathscr R=\frac43(1+2\arcsin(4/5)/\pi)
 +\frac2{\sqrt6}(1+2\arcsin(1/5)/\pi)
 +\frac1{\sqrt{21}}(1-2\arcsin(2/5)/\pi),
\]

\[
 \mathscr L=\frac4{3\sqrt3}
             +\frac1{\sqrt{15}}(1-2\arcsin(1/4)/\pi).
\]

Use $25/8<\pi<22/7$ and the rational bounds

\[
 \arcsin(4/5)<93/100,\quad \arcsin(1/5)<101/500,\quad
 \arcsin(2/5)>411/1000,\quad \arcsin(1/4)<253/1000.
\]

These follow respectively from the alternating sine bounds through
degrees $7,3,5,3$ on $[0,1]$. Squaring checks
$1/\sqrt3>577/1000$, $1/\sqrt{15}>258/1000$,
$1/\sqrt6<409/1000$, and $1/\sqrt{21}<219/1000$.
They imply

\[
 \mathscr R<\frac{2650234883}{825000000},\qquad
 \mathscr L>\frac{4619803}{4687500},\qquad
 G(1/2)<\frac{32058163928093}{60424804687500}<\frac9{16}.
 \tag{S36}
\]

For completeness the lower bound on $\pi$ follows already from
Machin's identity and
$\pi>16(1/5-(1/5)^3/3)-4/239>25/8$.
The classical upper bound follows from the positive integral
$\int_0^1x^4(1-x)^4/(1+x^2)dx=22/7-\pi$.
The checker verifies the displayed polynomial and rational inequalities;
the argument supplies their analytic connection.

**Full three-cell theorem in every ambient dimension.** Define

\[
 \kappa_3(\rho)=
 \begin{cases}
  \ell(1/\sqrt{1-\rho^2}-\sqrt3/2),&0<\rho\le1/2,\\
  \displaystyle\frac{\ell(1-\rho)(3\rho-1)}
                    {3\rho\sqrt{1-\rho^2}},&1/2<\rho<1.
 \end{cases}
\]

For every $d\ge2$, the cylindrical three-cell simplex and all balanced
Gaussian $L^2$ facet speeds satisfy

\[
 \boxed{\displaystyle
  \mathfrak q_\rho(f)\ge\kappa_3(\rho)
               \operatorname{dist}(f,\mathcal R_{3,d})^2,
                    \qquad0<\rho<1.}
 \tag{S37}
\]

Here $\mathcal R_{3,d}$ is the space of geometric rotation speeds;
$\kappa_3(\rho)>0$ is independent of $d$. In the active plane,
(S33)--(S35) prove the theorem. Their standard-sector bound is smaller
than $\ell\sigma$ on the applicable ranges, so it is valid for both
sectors. At $\rho=1/2$ the resulting full gap is at least
$\ell/(2\sqrt3)$.

To check every extra-dimensional mode, write
$x=(rp_{ij},z)$ with $z\in\mathbb R^{d-2}$ and decompose each speed
in the orthonormal Hermite basis in $z$. The surface measure factors as
$\nu\otimes\gamma_{d-2}$; the coefficient $a_\rho$ is unchanged and
the Mehler interaction in extra degree $m$ is multiplied by $\rho^m$.
Only degree $m=0$ has a mass constraint. In the active sign/standard
components the forms at degree $m$ are therefore

\[
 \rho^{-1}\left[\int a_\rho h^2d\nu
                       -\rho^m\iint K_\pm h h\,d\nu d\nu\right].
 \tag{S38}
\]

For sign components with $m\ge1$, (S18) implies
$K_-\le a_\rho$ as quadratic forms, so (S38) is at least
$\ell(1-\rho)\sigma^{-1}\|h\|^2$.
For standard degree $m=1$, (S38) is exactly the $K_+$ squared-difference
form, hence at least $\ell\sigma\inf_b\|h-b\|^2$ by (S31)--(S32).
Its kernel consists precisely of constant $h$.
For standard $m\ge2$, the positive-kernel row identity gives
$K_+\le a_\rho/\rho$ as forms, and again the gap is at least
$\ell(1-\rho)\sigma^{-1}$.

The only zero modes are thus the degree-zero sign mode $h=r$ and the
degree-one standard modes with constant $h$. They are exactly rotation
within the active plane and rotation mixing that plane with an unused
coordinate, respectively, of total dimension $1+2(d-2)=2d-3$.
Parseval and the bounded forms justify the orthogonal sum. The other
block gaps dominate $\kappa_3(\rho)$: in the upper range the ratio to
$\ell(1-\rho)/\sigma$ is $(3\rho-1)/(3\rho)<1$, and in the lower
range this follows from $\sigma\ge\sqrt3/2$. This proves (S37) for
all $L^2$ speeds, uniformly in $d$.

In particular, the continuation note's conditional measurable local
growth theorem is unconditional for three cells in the active plane at
every fixed positive correlation. If the initial three-cell global
optimality interval ends at $r_3<1$, its endpoint must have a nonsimplex
tied ordinary global maximizer: the additional-interface-zero-mode
alternative is now excluded. This remains a global rigidity obligation,
not a proof of full simplex optimality.

For $k\ge4$, (S4) supplies the precise symmetry blocks. The rotation hook
is settled below, while the other blocks still require estimates. The following reductions
separate the ambient-dimension issue, which can be settled for every $k$,
from the active-space spectral inequality, which remains open.

## Every cell count: the active gap lifts to all ambient dimensions

Fix $k\ge3$, $n=k-1$ and $0<\rho<1$. In this section the operators are
defined on the **unconstrained** active facet space $\mathcal H$ by the
same formulas as (B6). Write

\[
 L_\rho f=(a_{\rho,ij}/\rho)f_{ij},\qquad
 \langle f,\mathsf B_\rho f\rangle=B_\rho(f),\qquad
 \mu_\rho=\ell/\sigma,\qquad
 \mathcal T_k=\{t^v:t^v_{ij}=v\cdot n_{ij},\ v\in E\}.
 \tag{S39}
\]

The bounds and Hilbert--Schmidt proof in (B7) apply without a mass
constraint, so these are bounded self-adjoint operators and
$L_\rho\ge\mu_\rho I$. Equation (S5) gives
$\mathcal T_k=\operatorname{ran}C^*$ and
$\mathcal H=\mathcal T_k\mathbin\oplus\mathcal H_0$ orthogonally.

**Translation identity.** On the entire active space,

\[
       \mathsf B_\rho t=\rho^{-1}L_\rho t
       \quad(t\in\mathcal T_k),\qquad
       (L_\rho-\rho\mathsf B_\rho)t=0.
 \tag{S40}
\]

Here is a direct proof for the simplex, without a maximizing assumption.
The distributional gradient of $1_{T_i}$ is
$\sum_{j\ne i}n_{ij}\mathcal H^{n-1}|_{F_{ij}}$. Gaussian convolution
therefore gives
$\partial_vT_\rho1_{T_i}=\rho\sum_{j\ne i}\int
M_\rho(x,y)t^v_{ij}\,d\nu_{ij}(y)$.
Subtract the formula for $j$. Reflection swapping $i,j$ makes
$U_{\rho,ij}$ zero on the whole pair plane, so its gradient there is
$a_{\rho,ij}n_{ij}$. Comparing the two expressions yields
$\rho^2\mathsf B_\rho t=a_\rho t$, which is (S40).
All integrals converge by the Gaussian facet bounds. Equivalently one
can first use compact cutoffs in the distributional calculation and
then remove them. Heilman--Tarter, Lemma 5.4, contains the corresponding
translation formula for actual maximizers; their maximizing hypothesis
does not establish simplex stability. The direct calculation above
verifies the identity needed here.
[Primary source, Lemmas 4.6 and 5.4](https://stevenheilman.org/~heilman/papers/plurthreeB.pdf).

Assume only that the active balanced form is nonnegative. For any
$h=t+g$, $t\in\mathcal T_k$, $g\in\mathcal H_0$, (S40) gives

\[
\begin{split}
 \langle h,(L_\rho-\rho\mathsf B_\rho)h\rangle
 &=\rho\mathfrak q_\rho(g)+(1-\rho)\langle g,L_\rho g\rangle\\
 &\ge(1-\rho)\mu_\rho
                    \operatorname{dist}(h,\mathcal T_k)^2.
\end{split}
\tag{S41}
\]

In particular this unconstrained form is nonnegative with exactly the
translation kernel. The use of the orthogonal decomposition by $C^*$
is essential: no change from Gaussian mass orthogonality to
$a_\rho$-weighted orthogonality has been made.

**Dimension-lifting theorem.** Suppose, in the active space, that
$\mathfrak q_\rho(g)\ge\kappa_{\rm act}
\operatorname{dist}(g,\mathcal R_k)^2$ for all $g\in\mathcal H_0$,
where $\kappa_{\rm act}>0$. Then, for every $d\ge k-1$ and all balanced
Gaussian $L^2$ facet speeds of the cylindrical simplex,

\[
 \boxed{\displaystyle
 \mathfrak q_{\rho,d}(f)\ge
 \min\!\left\{\kappa_{\rm act},
              \frac{\ell(1-\rho)}{\sqrt{1-\rho^2}}\right\}
                  \operatorname{dist}(f,\mathcal R_{k,d})^2.}
 \tag{S42}
\]

To prove this, expand in Hermites in $z\in\mathbb R^{d-n}$.
In extra degree $m$ the active form is
$\mathfrak q_\rho^{[m]}(h)=\langle h,
(L_\rho-\rho^m\mathsf B_\rho)h\rangle$.
Only $m=0$ is mass constrained. Degree one is (S41). For $m\ge2$,

\[
 L_\rho-\rho^m\mathsf B_\rho
  =\rho^{m-1}(L_\rho-\rho\mathsf B_\rho)
                      +(1-\rho^{m-1})L_\rho
  \ge(1-\rho)\mu_\rho I.
 \tag{S43}
\]

The degree-zero kernel is $\mathcal R_k$. The degree-one kernel is
$\mathcal T_k$ tensored with the linear functions of $z$; these are
exactly rotations mixing the active and unused coordinates. Rotations
inside the unused space have zero facet speeds. Thus the full kernel
has dimension
$n(n-1)/2+n(d-n)=(k-1)(2d-k)/2$ and equals
$\mathcal R_{k,d}$. Orthogonal Hermite expansion and boundedness of the
forms justify the sum and prove (S42). For $k=2$, the active balanced
space is zero and the extra-degree forms are exactly
$\mu_\rho(1-\rho^{m-1})\|h\|_{\mathcal H}^2$: degree zero is forbidden, degree one
is rotation, and the sharp nonrotational gap for $d\ge2$ is
$\ell(1-\rho)/\sigma$.

Consequently the arbitrary-$k$ extension of (S37) is **equivalent to its
active-space version**. All extra-dimensional modes, their gaps and
their kernels are accounted for. In particular the known fixed-$k$
small-correlation active theorem already gives an ambient-dimension
independent gap on its interval. No active-space theorem at a new
correlation is inferred from this lifting argument.

## Every cell count: the reflected sign kernel is a determinant

Fix $k\ge3$, put $p=k-1$, and write

\[
 W=\{x\in E:x_1=x_2>x_3>\cdots>x_k\},\quad
 I=(2,3,\ldots,k),\quad \tau=\rho/\sigma^2,
\]

\[
 \Gamma_\rho(x,y)=\sigma^{-(k-1)}
  \exp\!\left[-\frac{\rho^2(|x|^2+|y|^2)}{2\sigma^2}\right].
\]

A sign-isotypic speed is determined by $h=f_{12}$, which is alternating
under permutations of $3,\ldots,k$. Restricting $h$ to $W$ makes it an
arbitrary scalar function; extend it alternatingly to $F_{12}$ and then
by the sign action to the other facets. Chamber walls have zero surface
measure. The full norm and quadratic form are both multiplied by
$\binom{k}{2}(k-2)!$ relative to this scalar description, so the same
factor cancels from spectral inequalities.

Let $\mathsf K_\rho=\rho\mathsf B_\rho$. Before antisymmetrizing the
remaining labels, its sign-sector kernel on $F_{12}$ is

\[
 \widetilde K_\rho(x,y)=2\Gamma_\rho(x,y)
 \left[e^{\tau x\cdot y}
        -\sum_{l=3}^k e^{\tau x\cdot(2l)y}\right].
 \tag{S44}
\]

Indeed edges disjoint from $12$ do not interact. For each $l$, the
$1l$ and $2l$ terms both have sign minus after applying the sign action
and their incidence signs. Their Mehler kernels agree on $H_{12}$,
because $(12)$ fixes $x$. The diagonal coefficient is $2$.
Also
$x\cdot(2l)y=x\cdot y-(x_2-x_l)(y_2-y_l)$.
For $k=3$ this is exactly $K_-$ in (S14).

For $k\ge4$ the kernel (S44) is negative near the origin:

\[
 \widetilde K_\rho(0,0)=2\sigma^{-(k-1)}(3-k)<0.
 \tag{S45}
\]

Continuity gives negativity on open relative-interior patches as well.
Thus the direct pointwise comparison of (S44) with a positive reflected
kernel fails. This is not a negative quadratic-form example: the
admissible functions on $F_{12}$ are alternating, and their cancellations
must be retained.

**Ordered-chamber formula.** After folding these cancellations, the exact
scalar kernel on $W$, integrating against $\nu_{12}|_W$, is

\[
 \boxed{\displaystyle
 \mathcal K_\rho(x,y)=
  2\Gamma_\rho(x,y)e^{\tau x_2y_2}
       \det\!\left[e^{\tau x_i y_j}\right]_{i,j\in I}>0
             \qquad(x,y\in W).}
 \tag{S46}
\]

For the identity, split the $y$ integral into chambers and sum (S44)
against the signs of $S_{\{3,\ldots,k\}}$. The representatives
$1,(23),\ldots,(2k)$, followed by this subgroup, enumerate every
permutation of $I$ once. Coordinate $1$ always contributes
$x_1y_1=x_2y_2$; the other coordinates give the determinant expansion,
with precisely its permutation signs. This proves (S46) for all $k$.

Here is a self-contained proof of its strict sign. For distinct ordered
real $a_i,b_j$ and $\tau>0$, the determinant
$\det[e^{\tau a_i b_j}]$ cannot vanish. Otherwise a nonzero combination
of $p$ exponentials with distinct real exponents would have $p$ distinct
zeros. Dividing by the exponential with least exponent, differentiating,
and using Rolle's theorem proves by induction that such a combination
has at most $p-1$ zeros. The determinant sign is constant on the connected
ordered domain. Expanding at $\tau=0$ determines this sign:

\[
 \det[e^{\tau a_i b_j}]
 =\frac{\tau^{p(p-1)/2}}{\prod_{j=0}^{p-1}j!}
    \prod_{i<j}(a_i-a_j)\prod_{i<j}(b_i-b_j)
                  +O(\tau^{p(p-1)/2+1})>0
 \tag{S47}
\]

for sufficiently small positive $\tau$, hence for every $\tau>0$.
This also identifies the same degree $d_k=p(p-1)/2$ as the Hermite
floor (S12). The finite checker compares the full coset sum with the
independent determinant expansion for $3\le k\le8$ and rejects omitted
or incorrectly signed adjacent terms. The universal proof is the coset
argument above, not this bounded replay.

The classical mechanism is the alternating method of images.
Grabiner, Theorem 5 and Section 4, equations (5)--(6), printed
pp. 184--186, give the reflection sum and Gaussian determinant for
Brownian motion killed on a Weyl-chamber wall. His definition on
pp. 182--183 requires the diffusion to be invariant under the reflection
group. Those hypotheses and the proof's path-cancellation step were
inspected. The determinant identity and positivity here were proved
directly; this comparison makes no claim of a new reflection principle.
[Primary source](https://www.numdam.org/article/AIHPB_1999__35_2_177_0.pdf).

## What the reflection mechanism still has to prove

For three cells the positive sign kernel has a known positive ground
state $r$, with $K_-r=a_\rho r$. This converts the form to the
nonnegative squared differences (S18). The standard sector instead has
the exact translation row identity, and the reflected even OU gap beats
its negative mass term in the overlapping correlation ranges. These
ground-state identities and quantitative comparisons, in addition to
pointwise positivity, are what prove (S37).

For $k\ge4$, (S46) supplies positivity in the sign chamber, but rotations
belong to the hook representation, not the sign representation. There
is no rotational sign ground state to substitute into (S18). The exact
sign-sector assertion is

\[
 \left\|a_\rho^{-1/2}\mathcal K_\rho
                         a_\rho^{-1/2}\right\|_{L^2(W,\nu_{12})}<1.
 \tag{S48}
\]

The operator is compact, self-adjoint and nonnegative as a quadratic
form, by restriction of the Mehler form. A positive kernel alone gives
neither this norm bound nor comparison with $a_\rho$. It is proved below
by (S61) for four cells and (S62)--(S66) for every higher cell count. A positive
supersolution $w$ satisfying
$\mathcal K_\rho w\le\theta_\rho a_\rho w$ with
$\theta_\rho<1$ suffices by the weighted Schur bound. The polynomial
supersolution in (S66) gives $\theta_\rho=\rho^{e_k-1}<1$ for $k\ge5$.

In particular (S46) cannot simply inherit the ordinary Weyl-chamber OU
spectrum. In the coordinates $a=(x_2,x_3,\ldots,x_k)$, the facet metric
and centering relation are
$2\,da_1^2+\sum_{j=2}^p da_j^2$ and
$2a_1+\sum_{j=2}^p a_j=0$. Permuting $a_1$ with another coordinate
preserves neither. The extra factor $e^{\tau x_2y_2}$ in (S46) records
this duplicated coordinate. Grabiner's reflection-invariance hypothesis
therefore does not identify this operator with the ordinary killed OU
semigroup in the facet metric. For $k=3$ the facet is one dimensional
and the separate pointwise comparisons (S32) supply the needed bridge.

The all-$k$ dimensional part of (S37) is settled by (S42). Its active
balanced gap for $k\ge4$, in the standard and remaining symmetry
blocks in (S4), is still an open mathematical obligation. The rotation
hook, two-row moment block and sign block are settled by
(S56), (S60) and (S62)--(S66) below. Neither a global simplex
inequality nor exclusion of tied maximizers follows from these
reductions.

## Every cell count: the standard sector is a scalar positive kernel

The standard block also has a scalar reduction for every $k\ge3$.
Formula (S4) associates $V_{(k-1,1)}$ only with the trivial
$S_{k-2}$ type. Thus its multiplicity functions $h$ on $F_{12}$ are
invariant under permutations of $3,\ldots,k$. Let $h_{ij}$ be the
transport of $h$ to $F_{ij}$; stabilizer invariance makes this definition
independent of the transporter. For $v\in E$ write
$f^{v,h}_{ij}=(v_i-v_j)h_{ij}$, using the same $h_{ji}=h_{ij}$.
These tensor speeds span the whole standard isotypic space. Their norms
and mass changes are

\[
 \|f^{v,h}\|^2=k|v|^2\|h\|_{L^2(\nu_{12})}^2,\qquad
       Cf^{v,h}=kv\int h\,d\nu_{12}.
\]

For a nonzero $v$, balance therefore imposes exactly the ordinary
$\nu_{12}$-mean-zero condition. The scalar interaction kernel on the
full reference facet is

\[
 K_{{\rm std},\rho}(x,y)=\Gamma_\rho(x,y)
   \left[2e^{\tau x\cdot y}
               +\sum_{l=3}^k e^{\tau x\cdot(2l)y}\right]>0,
 \qquad
 \int K_{{\rm std},\rho}(x,y)d\nu_{12}(y)=a_\rho(x)/\rho.
 \tag{S49}
\]

To verify the coefficients, pair the $1l$ and $2l$ terms in the $12$
row of (S10). Their geometric kernels agree after transport to the
reference facet, and their label coefficients subtract to
$(v_1-v_l)-(v_2-v_l)=v_1-v_2$. The diagonal coefficient is $2$.
This proves the kernel identity. Setting $h=1$ gives the row sum by
(S40), or directly by differentiating the reflected score. When a
particular $v_1-v_2$ is zero, use a different nonzero component or
linearity; the scalar operator does not depend on $v$.

Let $A=[e^{\tau x_i y_j}]_{i,j\in I}$ as in (S46), and let
$A_0=[e^{\tau x_i y_j}]_{i,j=3}^k$. Folding invariant functions to $W$
gives the exact positive kernel

\[
 \mathcal K_{{\rm std},\rho}(x,y)
  =\Gamma_\rho(x,y)e^{\tau x_2y_2}
       \left[\operatorname{per}A
                   +e^{\tau x_2y_2}\operatorname{per}A_0\right].
 \tag{S50}
\]

Indeed the first permanent enumerates every permutation of $I$ once,
while the second supplies one additional copy of the permutations fixing
$2$. This gives coefficient $2$ for the self term and $1$ for every
adjacent term, exactly as in (S49). No factorial normalization is
inserted: the kernel integrates against $\nu_{12}|_W$, and both scalar
forms and norms acquire the same factor $(k-2)!$ under unfolding.

Consequently the precise higher-cell extension of (S23) is the scalar
inequality below on invariant mean-zero $h$:

\[
\begin{split}
 \frac{\rho\mathfrak q_\rho(f^{v,h})}{k|v|^2}
 &=\int a_\rho h^2\,d\nu_{12}
               -\iint K_{{\rm std},\rho}(x,y)h(x)h(y)\,d\nu_{12}^2\\
 &=\frac12\iint K_{{\rm std},\rho}(x,y)(h(x)-h(y))^2d\nu_{12}^2
              -(1-\rho)\int\frac{a_\rho}{\rho}h^2d\nu_{12}.
\end{split}
\tag{S51}
\]

Strict positivity of (S51) on its stated constraint remains unproved
here for arbitrary correlation when $k\ge4$. Constants are the
translation ground state for the row-normalized kernel; they are
excluded by ordinary mass balance, whereas the squared-difference form
has the varying weight $a_\rho/\rho$. This is the same mass-versus-weight
issue identified after (S23). The permanent formula resolves the
reflection signs, not this weighted Poincare inequality. The hook is
treated next; the other blocks in (S4) require their own estimates as well.

## Every cell count: cone positivity settles the rotation hook

**Hook theorem.** For every $k\ge3$ and every $0<\rho<1$, the complete
$V_{(k-2,1,1)}$ isotypic part of the active interface form is
nonnegative. Its kernel is precisely the one geometric rotation copy,
and it has a positive $L^2$ gap modulo that copy. This includes mixing
between all multiplicity branches, and uses no maximizing hypothesis.
For $k=3$ this is the sign sector already treated above.

The useful generalization of a positive scalar ground state is a matrix
kernel preserving the ordered chamber. Its right ground state is the
rotation field $x$. A positive Gram matrix supplies a positive left
ground state as well. Both positivity statements are needed below.

### The complete hook kernel

Put $H=E\cap\{x_1=x_2\}$, $n=(e_1-e_2)/\sqrt2$ and
$G_0=S_{\{3,\ldots,k\}}$. The minus-eigenspace of $(12)$ in
$\Lambda^2E$ is exactly $n\wedge H$. Thus signed induction identifies
the hook multiplicity space with fields $v:F_{12}\to H$ satisfying
$v(\pi y)=\pi v(y)$ for $\pi\in G_0$, or, equivalently, arbitrary
$H$-valued fields on $W$. This is the whole multiplicity space, not
just fields of the form $h(x)x$. For example its scalar speeds are the
matrix coefficients of $g(n\wedge v(x))$ on the other facets, with a
fixed vector of $\Lambda^2E$ supplying the representation coordinate.
The inner product on the multiplicity space is the ordinary Euclidean
one on $H$, integrated against $\nu_{12}$; an overall positive constant
in the induced identification cancels from the form and the norm.

For $l\ge3$ define $w_l=(e_1+e_2)/2-e_l\in H$ and
$C_l=I_H-2w_lw_l^T$. The unnormalized interaction on these equivariant
fields is the matrix kernel

\[
 K_{{\rm hook},\rho}(x,y)=\Gamma_\rho(x,y)
 \left[2e^{\tau x\cdot y}I_H+
       \sum_{l=3}^k e^{\tau x\cdot(2l)y}C_l\right].
 \tag{S52}
\]

To check the coefficient, the adjacent transports $(2l)$ and $(1l)$
have equal Gaussian kernels on the reference plane. Both have incidence
coefficient $1$ when the transported edge keeps its orientation.
For $v\in H$, direct exterior algebra gives

\[
 (2l)(n\wedge v)+(1l)(n\wedge v)
                   =n\wedge\bigl(v-2w_l(w_l\cdot v)\bigr).
\]

Indeed put $a=e_2-e_l=w_l-n/\sqrt2$ and
$a'=e_1-e_l=w_l+n/\sqrt2$, and expand the two reflections.
The terms outside $n\wedge H$ cancel. This proves (S52), including
its coefficient $C_l$, whose eigenvalue along $w_l$ is $-2$.
Treating $C_l$ as an ordinary orthogonal reflection would be incorrect.

Let $r_m\in H$, $2\le m\le k-1$, have entries $k-m$ in the first
$m$ positions and $-m$ in the others. Write $V\eta=\sum_m\eta_mr_m$.
Then

\[
 (V^{-1}v)_j=(v_j-v_{j+1})/k,\quad 2\le j\le k-1,
 \qquad
 (V^TV)_{jm}=k\min(j,m)\bigl(k-\max(j,m)\bigr)>0.
\]

In particular $V$ carries the positive coordinate cone to $\overline W$.
After folding equivariant fields from $F_{12}$ to $W$, the kernel in
these coordinates is

\[
 \mathcal A_\rho(x,y)=V^{-1}
   \sum_{\pi\in G_0} K_{{\rm hook},\rho}(x,\pi y)\pi V.
 \tag{S53}
\]

### Positivity in chamber-gap coordinates

We prove the pointwise assertion

\[
       (\mathcal A_\rho(x,y))_{jm}>0
       \qquad(x,y\in W,\ 2\le j,m\le k-1).
 \tag{S54}
\]

The common factor $\Gamma_\rho$ is positive and can be suppressed.
For each permutation put $z=\pi y$, $v=\pi r_m$ and $b=k-m$.
Then $z_1=z_2=y_2$, $v_1=v_2=b$, and the remaining $v_i$ take two
values, $b$ and $-m$, in the same order as the corresponding $z_i$.
Write $E_0(z)=e^{\tau x\cdot z}$ and
$E_l(z)=e^{\tau x\cdot(2l)z}$.
The matrix entry is $1/k$ times the sum of the output differences
between coordinates $j$ and $j+1$ of
$2E_0v+\sum_lE_lC_lv$.

First let $j=2$. Put $d=b-v_3\ge0$. The self term together with the
$l=3$ term is
$2d(E_0(z)-E_3(z))\ge0$, because
$x\cdot z-x\cdot(23)z=(x_2-x_3)(z_2-z_3)>0$.
For each $l\ge4$, the remaining coefficient is $v_l-v_3$.
Pair $\pi$ with $(3l)\pi$. The coefficient changes sign, while

\[
 x\cdot(2l)z-x\cdot(2l)(3l)z
                         =(x_2-x_3)(z_l-z_3).
\]

The exponential difference therefore has the same sign as $v_l-v_3$;
each paired contribution is nonnegative. There is a permutation with
$v_3=-m$, since $m\le k-1$. Its self-minus-adjacent contribution is
strictly positive. This proves (S54) in row $2$.

Now let $j\ge3$ and put $d=v_j-v_{j+1}$. Pair each permutation with
$s\pi$, $s=(j,j+1)$. For the self term and each $l\notin\{j,j+1\}$,
the coefficient is respectively $2d$ and $d$. It changes sign under
the pairing, and the exponent difference has the sign of
$(x_j-x_{j+1})(z_j-z_{j+1})$. Thus all these contributions are
nonnegative. The two exceptional adjacent terms have coefficients
$D=2b-v_j-v_{j+1}\ge0$ and $-D$, respectively. Their paired sum is

\[
 D\bigl[E_j(z)-E_{j+1}(sz)+E_j(sz)-E_{j+1}(z)\bigr]\ge0,
\]

because the exponent differences in these two brackets are
$(x_j-x_{j+1})(z_2-z_{j+1})$ and
$(x_j-x_{j+1})(z_2-z_j)$, both strictly positive.
If $m\ge3$, the remaining coordinates include both values of $v_i$,
so the self pairing is strict for some permutation. If $m=2$, then
every remaining $v_i=-2$ and $D=2k>0$, making the exceptional pairing
strict instead. This proves every entry in (S54). Counting each pair
twice only inserts the same positive factor $1/2$ in these sums.

### The rotation field fixes the entire top eigenspace

Set $J=V^TV$ and work on $L^2(W,\nu_{12};\mathbb R^{k-2})$, first
with its equivalent $J$ inner product. The normalized operator is
$\mathcal T=a_\rho^{-1/2}\mathcal A_\rho a_\rho^{-1/2}$.
It is compact and self-adjoint in this inner product, by (B7) and the
original symmetric interface form. In ordinary coordinate inner
product, $\mathcal T^*J=J\mathcal T$.

Rotational equivariance gives $K_{{\rm hook},\rho}x=a_\rho x$.
For completeness, differentiating the rotated indicators gives
$\delta 1_{T_i}=\lambda_i^{f^\Omega}$, whereas commutation of OU with
rotations gives
$\delta U_{\rho,ij}=-(\Omega x)\cdot\nabla U_{\rho,ij}
=a_{\rho,ij}f^\Omega_{ij}$ on the facet. Thus
$\mathsf K_\rho f^\Omega=a_\rho f^\Omega$ without a mass
restriction. In hook multiplicity coordinates this is the claimed
identity for $x$.
Consequently

\[
 \psi(x)=\sqrt{a_\rho(x)}V^{-1}x>0,\qquad
 \chi(x)=J\psi(x)>0,\qquad
       \mathcal T\psi=\psi,\quad\mathcal T^*\chi=\chi.
 \tag{S55}
\]

Both fields are in $L^2$, and every entry of $J$ is positive. The
left ground state cannot be omitted: the gap-coordinate basis is not
orthogonal, so entrywise positivity is not itself self-adjointness in
the ordinary coordinate norm.

Here is the required positive-operator argument directly. For every
nonzero $L^2$ eigenvector $u$ with eigenvalue $\lambda$, positivity of
the matrix kernel and Tonelli's theorem give

\[
 |\lambda|\int\chi\cdot|u|\,d\nu_{12}
 \le\int\chi\cdot\mathcal T|u|\,d\nu_{12}
 =\int\chi\cdot|u|\,d\nu_{12}.
\]

The common integral is finite and strictly positive by Cauchy--Schwarz
and (S55). Hence $|\lambda|\le1$. Compactness accounts for every
nonzero spectral value. If a real eigenvector has eigenvalue $1$,
equality in the triangle inequality, together with the strictly
positive entries (S54), forces all its components to have one common
sign almost everywhere. Two independent eigenvectors at $1$ would
give a nonzero one orthogonal to $\psi$ in the $J$ inner product.
That is impossible for a one-sign field, since $J\psi>0$. The
eigenvalue $1$ is therefore simple.

Let $\delta_{k,\rho}=1-\sup\bigl(\{0\}\cup
\operatorname{spec}(\mathcal T|_{\psi^{\perp_J}})\bigr)>0$.
Conjugating back and using $a_\rho/\rho\ge\ell/\sigma$ proves

\[
 \boxed{\displaystyle
 \mathfrak q_\rho(f)\ge\frac{\ell}{\sigma}\delta_{k,\rho}
             \operatorname{dist}(f,\mathcal R_k)^2
       \quad\text{on the entire rotation hook}.}
 \tag{S56}
\]

The multiplicity ground state $v(x)=x$ supplies exactly one
$\Lambda^2E$ copy in the full isotypic space, namely geometric
rotations. This proves the hook theorem. The constant is a proved
positive spectral constant, not an evaluated explicit function of
$k,\rho$. The argument is a matrix version of the positive
ground-state mechanism behind (S18). The two-row and sign blocks are
treated next; the standard and remaining higher-cell
remaining block bounds require further arguments.

## Every cell count: a strict supersolution for the two-row moment sector

**Two-row theorem.** For every $k\ge4$ and every $0<\rho<1$, the
complete $V_{(k-2,2)}$ isotypic part has a strictly positive $L^2$ gap.
It has no zero modes. The argument uses the hook identity, but does not
assume positivity in the standard sector.

Put $q=k-2$ and
$H^0=\{v\in E:v_1=v_2=0,\ \sum_{l=3}^kv_l=0\}$.
In $\operatorname{Sym}^2E=\mathbf1\oplus E\oplus V_{(k-2,2)}$,
the $(12)$ minus-space is $n\odot H$. Its standard summand is the
line $n\odot p$, where $p=(e_1+e_2)/2-\mathbf1/k$, and its
orthogonal complement is $n\odot H^0$. Here $u\odot v$ denotes
$uv^T+vu^T$. To check the standard line, the equivariant map
$z\mapsto P_E\operatorname{diag}(z)P_E$ sends $n$ to $n\odot p$.
Thus the entire two-row multiplicity space consists of $H^0$-valued
$G_0$-equivariant fields on $F_{12}$, or arbitrary fields on $W$.

The symmetric, rather than exterior, transport identity is

\[
 (2l)(n\odot v)+(1l)(n\odot v)=n\odot v.
\]

It follows by expanding the same two reflections as in (S52). The
quadratic $a\odot a$ terms now contribute and cancel the
$-2n\odot w_l(w_l\cdot v)$ term. Consequently the full-reference
kernel in this sector is simply

\[
             K_{{\rm std},\rho}(x,y)I_{H^0}.
 \tag{S57}
\]

For $1\le m\le q-1$, let $s_m\in H^0$ have first two entries zero,
the next $m$ entries $q-m$, and the rest $-m$. Let $V_0$ have these
columns and $J_0=V_0^TV_0$; its entries are
$q\min(j,m)(q-\max(j,m))>0$.
The folded kernel

\[
 \mathcal A^{(2)}_\rho(x,y)=V_0^{-1}
       \sum_{\pi\in G_0}K_{{\rm std},\rho}(x,\pi y)\pi V_0
                   \quad\text{is strictly entrywise positive.}
 \tag{S58}
\]

Here is the additional pairing needed to prove this assertion. As before,
put $z=\pi y$, $v=\pi s_m$, and consider an output difference between
coordinates $j,j+1$, $j\ge3$. Every coefficient in (S57) is
$d=v_j-v_{j+1}$. Pair with $s\pi$, $s=(j,j+1)$.
The self term and the terms $l\notin\{j,j+1\}$ are nonnegative by
the same ordered-exponential comparison as in (S54). For the two
exceptional terms the paired contribution is
$d\Phi$, where
$\Phi=E_j(z)+E_{j+1}(z)-E_j(sz)-E_{j+1}(sz)$.
Choose the member of the pair with $z_j>z_{j+1}$, so $d\ge0$.
Writing

\[
 \alpha=\tau(x_2-x_j)\ge0,\quad
 \beta=\tau(x_j-x_{j+1})>0,\quad
 u=z_2-z_j>0,\quad t=z_j-z_{j+1}>0,
\]

and dividing by $E_j(sz)>0$, direct expansion gives

\[
 \frac{\Phi}{E_j(sz)}
 =(e^{\alpha t}-1)(e^{\beta t}-e^{-\beta u})
                          +(e^{\beta t}-1)>0.
\]

The self term is strict for some permutation for every input ray
$s_m$, because both values of its remaining coordinates occur.
This proves (S58) for all $k$.

Let $P_0$ be orthogonal projection onto $H^0$ and put
$v_0(x)=P_0x$. Projecting the rotation identity in (S52) onto $H^0$
gives the exact strict-supersolution residual

\[
 a_\rho v_0-K_{{\rm std},\rho}v_0
       =2\sum_{l=3}^k\left(e_l-\frac1q\sum_{r=3}^ke_r\right)R_l(x),
 \qquad
 R_l(x)=\int_{F_{12}}\Gamma_\rho(x,y)
    e^{\tau x\cdot(2l)y}(y_2-y_l)d\nu_{12}(y).
 \tag{S59}
\]

For $j\ge3$, substitute $(j,j+1)y$ into $R_{j+1}$. The difference
of the two exponents is
$\tau(x_j-x_{j+1})(y_2-y_{j+1})>0$ in the facet interior.
Both integrals then have the same nonnegative factor $y_2-y_j$.
It follows that $R_j(x)>R_{j+1}(x)$ for every $x\in W$.
Thus every gap-coordinate component of (S59) is strictly positive.

Normalize $\mathcal A^{(2)}$ by $a_\rho^{-1/2}$ as before, obtaining
$\mathcal T_2$, self-adjoint in the $J_0$ inner product. The field
$\psi_2=\sqrt{a_\rho}V_0^{-1}v_0$ is strictly positive, and (S59)
says $\mathcal T_2\psi_2<\psi_2$ componentwise. Therefore
$\chi_2=J_0\psi_2>0$ satisfies
$\mathcal T_2^*\chi_2<\chi_2$, since every entry of $J_0$ is positive.
All these fields are in $L^2$, by the Gaussian first-moment bounds.
For a nonzero eigenvector $u$ with eigenvalue $\lambda$, the direct
argument used after (S55) now gives the strict inequality

\[
 |\lambda|\int\chi_2\cdot|u|d\nu_{12}
 \le\int\mathcal T_2^*\chi_2\cdot|u|d\nu_{12}
 <\int\chi_2\cdot|u|d\nu_{12}.
\]

Hence $|\lambda|<1$ for every nonzero spectral value. Compactness
gives $\delta^{(2)}_{k,\rho}=1-\sup\operatorname{spec}\mathcal T_2>0$,
and

\[
 \boxed{\displaystyle
 \mathfrak q_\rho(f)\ge\frac{\ell}{\sigma}
       \delta^{(2)}_{k,\rho}\|f\|^2
              \quad(f\text{ in }V_{(k-2,2)}).}
 \tag{S60}
\]

These are positive spectral constants; their numerical values have not
been evaluated. This theorem treats the entire two-row block, not only
the finite-dimensional part seen by the moment map.

## Four cells: domination settles the remaining sign block

For $k=4$, $H^0$ is a line and (S57) acts on functions odd under
$(34)$. Let $S=\Gamma_\rho\sum_{\pi\in G_0}\operatorname{sgn}(\pi)
e^{\tau x\cdot\pi y}$ and
$U=\Gamma_\rho\sum_{\pi\in G_0}\operatorname{sgn}(\pi)
\sum_{l=3}^4e^{\tau x\cdot(2l)\pi y}$.
The two-row and sign scalar kernels on $W$ are respectively $2S+U$
and $2S-2U$. The exceptional pairing just proved, with $j=3$,
shows $U>0$ for $x,y\in W$. Consequently

\[
        0<\mathcal K_{{\rm sign},\rho}
                  <\mathcal K_{(2,2),\rho}\quad\text{on }W\times W.
 \tag{S61}
\]

The positive function $w(x)=x_3-x_4$ is a strict supersolution for
the larger kernel by (S59), and therefore for the sign kernel too.
After normalization by $a_\rho^{-1/2}$ it is a positive strict left
supersolution as well, since this scalar kernel is symmetric.
The same eigenvalue argument proves (S48) for $k=4$: its sign sector
has no zero modes and has a strictly positive $L^2$ gap for every
$0<\rho<1$. The higher-cell sign sectors require the separate
Gaussian-matrix comparison (S62)--(S66).

Combining (S56), (S60) and (S61), the **entire four-cell active-space
interface assertion is equivalent to the single standard-sector
inequality (S51)**. A proof of that inequality, with a positive gap,
would give the full four-cell analogue of (S37) in every ambient
dimension by (S42). That scalar inequality is not proved here.

## Every cell count: a Vandermonde supersolution settles the sign sector

**All-cell sign theorem.** For every $k\ge4$ and every $0<\rho<1$,
the entire active sign sector has a strictly positive $L^2$ gap and
no kernel. For $k\ge5$ the following bound is explicit. Set
$q=k-2$ and $e_k=q(q-1)/2$. Then

\[
 \mathfrak q_\rho(f)\ge\frac{\ell}{\sigma}
                  (1-\rho^{e_k-1})\|f\|^2
       \qquad(k\ge5,\ f\text{ in the sign sector}).
 \tag{S62}
\]

For $k=4$ the positive spectral gap from (S61) applies. The proof below
also supplies a strict supersolution in that case. Three-cell sign
speeds instead have their geometric rotation kernel, as in (S18).

### A totally positive minor comparison

If $A$ is a strictly totally positive $p\times p$ matrix, $p\ge2$,
and $D=A_{2:p,2:p}$, then

\[
                   0<\det A<A_{11}\det D.
 \tag{S63}
\]

Here is an elementary proof, including the needed hypothesis. Eliminate
the last row and column, using its positive diagonal entry as pivot.
Every minor of the resulting Schur complement is a bordered minor of
$A$ divided by that pivot, so the complement is strictly totally
positive. Its top-left entry is strictly smaller than $A_{11}$.
The ratio $\det A/\det D$ is unchanged by this elimination. Repeat
until the ratio is a scalar pivot; this proves (S63). Positive entries
alone would not suffice for this inequality.

Apply (S63) to $A=[e^{\tau x_i y_j}]_{i,j=2}^k$. Every one of its
minors is positive by the exponential argument used in (S47).
Writing
$D=[e^{\tau x_i y_j}]_{i,j=3}^k$, define the folded self-facet kernel

\[
 S_\rho(x,y)=\Gamma_\rho(x,y)e^{2\tau x_2y_2}\det D.
 \qquad
        0<\mathcal K_\rho(x,y)<2S_\rho(x,y)
                    \quad(x,y\in W).
 \tag{S64}
\]

This $S_\rho$ is precisely the alternating fold, in labels $3,\ldots,k$,
of the single same-facet Mehler kernel. No adjacent-facet terms or
factorials are included in it.

### The Gaussian matrix comparison

Put $w(x)=\Delta(x_3,\ldots,x_k)
=\prod_{3\le i<j\le k}(x_i-x_j)>0$ on $W$, and let
$R_{\rm self}(x)=\int_{F_{12}}M_\rho(x,y)d\nu_{12}(y)$.
We will prove

\[
 \int_W S_\rho(x,y)w(y)d\nu_{12}(y)
             \le\rho^{e_k}w(x)R_{\rm self}(x).
 \tag{S65}
\]

First, the restricted Mehler measure equals $\phi(0)/\sigma$ times
the law of a Gaussian on $H$ with mean $\rho x$ and covariance
$\sigma^2I_H$. For functions invariant under common translation of
all coordinates, that law can be represented by independent variables

\[
 V\sim N(\rho x_2,\sigma^2/2),\qquad
 Z_l\sim N(\rho x_l,\sigma^2),\quad 3\le l\le k.
\]

Indeed subtract $(2V+\sum_lZ_l)/k$ from all $k$ coordinates of
$(V,V,Z_3,\ldots,Z_k)$. The resulting vector has the stated Gaussian
law on $H$. Both the event $V\ge\max_lZ_l$ and the Vandermonde of
the $Z_l$ are invariant under this subtraction. Hence

\[
 R_{\rm self}=\frac{\phi(0)}\sigma
                 \mathbb P\{\max_lZ_l\le V\},\qquad
 \int_WS_\rho w\,d\nu_{12}
  =\frac{\phi(0)}\sigma
          \mathbb E\bigl[\Delta(Z)\mathbf1_{\{\max_lZ_l\le V\}}\bigr].
\]

The second identity follows either by splitting the full facet into
remaining-coordinate chambers or directly from the determinant fold.
The Vandermonde is signed on the full facet, so it cannot be bounded
pointwise by its value at the mean.

The needed imported identity is the Hermitian-matrix realization of
the Vandermonde transform. For ordered $b=(x_3,\ldots,x_k)$, let
$G$ be a full $q\times q$ Hermitian Gaussian matrix with independent
real diagonal entries of variance $1$, and independent real and
imaginary parts of each upper off-diagonal entry of variance $1/2$.
Then the ordered eigenvalues of
$\mathsf H=\rho\operatorname{diag}(b)+\sigma G$ have density

\[
 \det[\phi_\sigma(y_i-\rho b_j)]_{i,j=1}^q
                \frac{\Delta(y)}{\rho^{e_k}\Delta(b)}
       \quad\text{on }y_1>\cdots>y_q.
\]

This is Grabiner, Theorem 3, combined with the determinant formula (6),
Section 6.1 and the Vandermonde $h$-transform in Section 6.2,
printed pp. 179, 186 and 197--200. The Gaussian normalization on
p. 197, the arbitrary-start $h$-transform on p. 198, and harmonicity
in Theorem 7 on p. 199 were inspected in the
[retained primary text](sources/grabiner-weyl-chamber-1999.pdf).
Apply the Brownian theorem at time $\sigma^2$, starting at the
ordered vector $\rho b$. Its denominator is
$\Delta(\rho b)=\rho^{e_k}\Delta(b)$. We use the full Hermitian
space, including its trace coordinate; replacing it by the traceless
space would change the diagonal covariance. This is an application
of the classical matrix theorem, not a new noncollision theorem.

Split the independent Gaussian integral into ordered chambers and
antisymmetrize. The displayed density gives, for every real $v$,

\[
 \mathbb E[\Delta(Z)\mathbf1_{\{\max Z_l\le v\}}]
    =\rho^{e_k}\Delta(b)\,
                   \mathbb P\{\lambda_{\max}(\mathsf H)\le v\}.
\]

The diagonal of $\mathsf H$ has exactly the independent law of the
$Z_l$. For every Hermitian matrix its largest eigenvalue is at least
each diagonal entry, by the Rayleigh quotient on the coordinate
vectors. Couple $\mathsf H$ with its own diagonal and with the
independent $V$ above. Integrating the last display against $V$ proves
(S65). All expectations are finite Gaussian polynomial moments.

### The spectral bound

The row identity (S49) gives
$a_\rho=\rho(2R_{\rm self}+\sum_{l=3}^kR_l^{\rm adj})$, where
every adjacent row integral $R_l^{\rm adj}$ is strictly positive.
Combining this with (S64)--(S65) gives the pointwise supersolution

\[
         0<\mathcal K_\rho w
              <2\rho^{e_k}R_{\rm self}w
              <\rho^{e_k-1}a_\rho w
                     \quad\text{on }W.
 \tag{S66}
\]

For $k\ge5$, $\theta=\rho^{e_k-1}<1$. The compact symmetric
normalized sign operator $T=a_\rho^{-1/2}\mathcal K_\rho
a_\rho^{-1/2}$ has the positive $L^2$ supersolution
$\psi=\sqrt{a_\rho}w$ with $T\psi\le\theta\psi$.
For an eigenvector $u$, integrate $|Tu|\le T|u|$ against $\psi$,
as after (S55), to obtain $|\lambda|\le\theta$.
Conjugating back and using $a_\rho/\rho\ge\ell/\sigma$ proves (S62).
For $k=4$, $e_k=1$ and the same argument with the strict inequality
in (S66) gives $|\lambda|<1$ for every nonzero eigenvalue;
compactness then supplies a positive gap.

This comparison uses only the ordinary Weyl-chamber motion of the
untied $q$ coordinates. The duplicated coordinate remains the
independent scalar threshold $V$ with variance $\sigma^2/2$. It does
not identify the whole facet operator with an ordinary Weyl-chamber
OU semigroup, so the metric obstruction noted after (S48) is respected.

## Every symmetry block is a compressed totally positive immanant kernel

There is one formula behind all of the scalar and matrix kernels above.
For a real irreducible orthogonal representation $R_\lambda$ of $S_k$,
put $s=(12)$, $P_\lambda=(I-R_\lambda(s))/2$ and
$V_\lambda^-=\operatorname{ran}P_\lambda$.
The complete $\lambda$ multiplicity space is
$L^2(W,\nu_{12};V_\lambda^-)$: ordering all untied coordinates leaves
only the oriented-facet stabilizer $\langle s\rangle$.
Schur averaging gives its ordinary vector inner product, up to the
constant $k!/(2\dim V_\lambda)$, which cancels from spectral bounds.

For $A=[e^{\tau x_i y_j}]_{i,j=2}^k$ define an element of
$\mathbb R[S_{\{2,\ldots,k\}}]$ by

\[
 Z(A)=\sum_{\pi\in S_{\{2,\ldots,k\}}}
                   \left(\prod_{i=2}^kA_{i,\pi(i)}\right)\pi^{-1}.
\]

Then the full unnormalized interaction kernel in this block is exactly

\[
 \boxed{\displaystyle
 \mathcal K_{\lambda,\rho}(x,y)=
 2\Gamma_\rho(x,y)e^{\tau x_2y_2}
       \left.P_\lambda R_\lambda(Z(A))P_\lambda
                                      \right|_{V_\lambda^-}.}
 \tag{S67}
\]

In particular the block need not be scalar, but its finite matrix size
is exactly $\dim V_\lambda^-$; its entries are explicit immanants of
the same strictly totally positive exponential matrix.

Here is the group-algebra proof, including the normalization. Put
$d=e_1-e_2$, $m_g=e^{\tau x\cdot gy}$, and
$W_1=\sum_{g:g(1)=1}m_g g=e^{\tau x_2y_2}Z(A)$.
Because $x_1=x_2$ and $y_1=y_2$, $m_g$ is unchanged under left or
right multiplication by $s$. Thus direct coefficient comparison gives

\[
 \sum_{g\in S_k}(d\cdot gd)m_g g=(1-s)W_1(1-s).
 \tag{S68}
\]

For example the coefficient of $g$ on either side is $m_g$ times
$\mathbf1_{g(1)=1}-\mathbf1_{g(1)=2}
-\mathbf1_{g(2)=1}+\mathbf1_{g(2)=2}$.
The signed facet kernel sums over cosets of $\langle s\rangle$;
the two representatives of each coset contribute equally after acting
on $V_\lambda^-$. Therefore the full-group sum is divided by $2$.
Each factor $1-s$ then acts as $2P_\lambda$, giving the factor $2$
in (S67), rather than $1$ or $4$.
The sign representation reduces (S67) to (S46). The standard
representation reduces it to (S50), including its additional fixed-row
permanent. The exterior and symmetric second powers give (S52) and
(S57).

**Pointwise positivity in every block.** There is a fixed orthonormal
Young basis, independent of $x,y,\rho$, in which every entry of
$\mathcal K_{\lambda,\rho}(x,y)$ is nonnegative for $x,y\in W$.
This is an entrywise assertion, not the interface stability inequality.
Two classical representation and total-positivity inputs suffice:

* In Young's orthogonal form, each adjacent transposition has scalar
  blocks $1,-1$ or a two-dimensional block with diagonal entries
  $1/r,-1/r$ and both off-diagonal entries
  $\sqrt{1-r^{-2}}$, where $|r|\ge2$. Consequently $I+R_\lambda(s_i)$
  is entrywise nonnegative. Also $R_\lambda((12))$ is diagonal: the
  entries $1,2$ of a standard tableau share a row or a column.
  See Vershik--Okounkov, Propositions 6.1--6.2 and (6.3)--(6.5),
  [retained arXiv v3, p. 23](sources/vershik-okounkov-young-form-0503040v3.pdf).
* For every totally nonnegative $A$, the group-algebra sum
  $\sum_\pi(\prod_iA_{i,\pi(i)})\pi$ is a nonnegative linear
  combination of products of interval symmetrizers
  $z_J=\sum_{\pi\in S_J}\pi$, with $J$ consecutive.
  Rhoades--Skandera state this exact decomposition in the proof of
  Proposition 2.2 and describe its cone inclusion in Section 7,
  [retained author preprint, pp. 4 and 15](sources/rhoades-skandera-kl-immanants-2005.pdf),
  crediting Stembridge's Theorem 2.1. Inverting group elements reverses
  the products, and each $z_J$ is invariant under inversion, so the
  same statement holds for the convention defining $Z(A)$ above.

The cited statements and relevant formulas were inspected in the
retained primary texts. Their hypotheses match because every minor of
$A$ is strictly positive by (S47). No positivity result for arbitrary
positive-entry matrices is imported.

For completeness, the passage from these inputs to matrix positivity
does not require a further representation theorem. On an interval
$J$, average the orthogonal projections
$(I+R_\lambda(s_i))/2$ over its adjacent generators. This is an
entrywise nonnegative self-adjoint operator with spectrum in $[0,1]$.
Its powers converge to the orthogonal projection onto the common
fixed space, which is $R_\lambda(z_J)/|S_J|$. Hence every interval
symmetrizer, and every product of them, is entrywise nonnegative in
the same Young basis. The decomposition above proves this for
$R_\lambda(Z(A))$. Finally $P_\lambda$ is a diagonal $0,1$ matrix in
that basis, so compression preserves entrywise nonnegativity. Thus

\[
        (\mathcal K_{\lambda,\rho}(x,y))_{ab}\ge0
              \qquad(x,y\in W)
 \tag{S69}
\]

for every nontrivial $\lambda$. No assertion of strict positivity of
all entries is needed here.

This provides an all-$k$ structural explanation of the reflected
kernel signs. It does not place the normalized spectrum below $1$.
Indeed the standard block already has the positive translation
eigenfunction with eigenvalue $1/\rho>1$ before imposing mass balance.
The remaining tasks are still the constrained standard inequality
(S51), and appropriate spectral bounds or positive supersolutions for
the other blocks when $k\ge5$. Equations (S56), (S60) and
(S62) provide those additional arguments for the sectors they cover.

## The exact remaining standard-sector spectral test

There are two conditions to check in the standard block, even after
its positive translation eigenfunction is known. This distinction is
essential when using the reflected kernels to seek an all-$k$ version
of (S37).

Work on $L^2(W,\nu_{12})$ with the folded kernel (S50), and put

\[
 T=a_\rho^{-1/2}\mathcal K_{{\rm std},\rho}a_\rho^{-1/2},\qquad
 \psi=\sqrt{a_\rho},\quad b=a_\rho^{-1/2},\quad
 e=\psi/\|\psi\|,\quad c=\rho^{-1}-1>0.
\]

The row identity gives $T\psi=\rho^{-1}\psi$. Since the kernel is
strictly positive, this is the simple largest eigenvalue. One can check
the assertion directly: for an eigenvector $v$, positivity gives
$|\lambda|\langle\psi,|v|\rangle\le
\langle T\psi,|v|\rangle=\rho^{-1}\langle\psi,|v|\rangle$;
equality at the positive eigenvalue forces one sign, and orthogonality
to $\psi$ then rules out another eigenvector. Compactness supplies the
remaining spectral separation. Define

\[
 \lambda_* =\sup\operatorname{spec}(T|_{e^\perp}),\qquad
 b_1=\langle b,e\rangle>0,\quad b_\perp=b-b_1e.
\]

For $\lambda_*<1$, let $A=(I-T)|_{e^\perp}>0$. The standard
inequality (S51), with a positive gap, is equivalent to the two tests

\[
 \boxed{\quad \lambda_*<1,\qquad
       c\langle b_\perp,A^{-1}b_\perp\rangle<b_1^2.\quad}
 \tag{S70}
\]

Indeed $u=\sqrt{a_\rho}h$ converts ordinary mass balance into
$\langle u,b\rangle=0$. Write $u=te+v$ with $v\perp e$.
Then $t=-\langle v,b_\perp\rangle/b_1$ and

\[
 \langle u,(I-T)u\rangle
   =\langle v,Av\rangle
        -\frac{c}{b_1^2}\langle v,b_\perp\rangle^2.
\]

If $\lambda_*\ge1$, compactness supplies an eigenvector at an
eigenvalue at least $1$ in $e^\perp$. Its span with $e$ intersects
$b^\perp$ nontrivially, and the form is nonpositive there. Thus a
positive balanced gap requires $\lambda_*<1$. Under that condition,
Cauchy--Schwarz in the $A$ inner product proves sufficiency of the
second test. It also proves necessity, by taking
$v=A^{-1}b_\perp$ if $b_\perp\ne0$. Equality in the second test
gives a balanced zero mode, and reversing it gives a negative mode.
The case $b_\perp=0$ is immediate.

This also gives a quantitative implication. Put
$\beta=c\langle b_\perp,A^{-1}b_\perp\rangle/b_1^2$ and
$\delta=1-\lambda_*>0$. The normalized form on $b^\perp$ is at
least

\[
 (1-\beta)\delta\frac{b_1^2}{\|b\|^2}\|u\|^2.
\]

Consequently the standard part of $\mathfrak q_\rho$ has gap at least
$\mu_\rho(1-\beta)\delta b_1^2/\|b\|^2$ if (S70) is proved.
Equivalently its second test is
$\langle b,(I-T)^{-1}b\rangle<0$, after the first test makes this
inverse well defined. Neither test has been established here for the
simplex at every $k\ge4$ and every $0<\rho<1$.

**A control for the mass-weight mismatch.** On two atoms of measure
$1/2$ each take $\rho=1/2$, $a=(1,4)$ and the symmetric positive
kernel matrix

\[
 K=\begin{pmatrix}2&2\\2&14\end{pmatrix}.
\]

Its integral row sum is $a/\rho$. The normalized integral operator
has matrix
$T=\left(\begin{smallmatrix}1&1/2\\1/2&7/4\end{smallmatrix}\right)$,
eigenvalues $2$ and $3/4$, and positive translation eigenvector
$\sqrt a=(1,2)$. Nevertheless the ordinarily balanced $h=(1,-1)$
satisfies

\[
 \int ah^2d\nu-\iint K(x,y)h(x)h(y)d\nu(x)d\nu(y)=-\tfrac12.
\]

Here $b=(1,1/2)$ and
$\langle b,(I-T)^{-1}b\rangle=1/2>0$; the second test fails.
This is a counterexample to discarding the scalar test or replacing
ordinary mass balance by spectral orthogonality. It is not a simplex
counterexample. The checker verifies these two-atom identities exactly.

## The standard and rotation-hook kernels have an exact Gaussian link

There is also a differential explanation of the hook transport matrix
and the rotational ground state. It keeps the junction term visible
when trying to extend the one-dimensional OU comparison.

Use the full-facet scalar kernel (S49) and vector hook kernel (S52),
denoted $K_{\rm std}$ and $\mathsf K_{\rm hook}$. Fix $k\ge3$ and
$0<\rho<1$. Let $h$ be a $C^1$ function on
$H=E\cap\{y_1=y_2\}$, invariant under $S_{\{3,\ldots,k\}}$,
such that $h$ and its first derivatives have at most polynomial
growth. Gradients below are intrinsic to $H$. Put

\[
 w_l=(e_1+e_2)/2-e_l,\quad D_l=I-w_lw_l^T,\quad
 C_l=I-2w_lw_l^T,\qquad
 (\mathcal R_lh)(x)=\int_{F_{12}}\Gamma_\rho(x,y)
    e^{\tau x\cdot(2l)y}(y_2-y_l)h(y)d\nu_{12}(y).
\]

Here $D_l=P_H(2l)|_H$ is the tangent compression of the adjacent
transposition, whereas $C_l$ is its hook transport. Since
$|w_l|^2=3/2$, they satisfy $C_lD_l=D_lC_l=I$ on $H$.
In particular $C_l$ is the inverse tangent compression, rather than
an orthogonal reflection. The following identities hold on the
interior of $F_{12}$:

\[
 \boxed{\begin{aligned}
 \mathsf K_{\rm hook}(y h(y))(x)
   &=\sigma^2\mathsf K_{\rm hook}(\nabla h)(x)
                           +\rho x(K_{\rm std}h)(x),\\
 \nabla(K_{\rm std}h)(x)
   &=\rho\mathsf K_{\rm hook}(\nabla h)(x)
                   +\frac{\rho}{\sigma^2}
                                    \sum_{l=3}^k w_l\mathcal R_lh(x).
 \end{aligned}}
 \tag{S71}
\]

To prove the first identity, integrate each term of
$\mathsf K_{\rm hook}\nabla h$ by parts in $y$. Write
$M_0=\Gamma_\rho e^{\tau x\cdot y}$ and
$M_l=\Gamma_\rho e^{\tau x\cdot(2l)y}$. Gaussian differentiation
gives

\[
 \nabla_y(M_0\phi_H)=-\sigma^{-2}(y-\rho x)M_0\phi_H,
 \qquad
 C_l\nabla_y(M_l\phi_H)
            =-\sigma^{-2}(C_ly-\rho x)M_l\phi_H.
\]

The constant factor between $\phi_Hd\mathcal H^{k-2}$ and
$d\nu_{12}$ cancels. There are no surviving boundary terms.
Indeed on the face $y_2=y_l$ the outward unit normal is
$m_l=-w_l/|w_l|$, and $M_l=M_0$. Thus the same-facet term and
the $l$th adjacent term cancel there, since
$(2I+C_l)m_l=0$. For $r\ne l$,
$C_rw_l=w_l-w_r=-C_lw_r$. Substitution $y\mapsto(lr)y$
maps the face $l$ to face $r$, preserves $h$ and its Gaussian surface
measure, and converts $M_r$ into $M_l$: on the face $l$,
$(2r)y=(2l)(lr)y$. These terms cancel in pairs. Gaussian decay
removes the boundary at infinity. The remaining interior integral is

\[
 \mathsf K_{\rm hook}\nabla h
    =\sigma^{-2}\bigl(\mathsf K_{\rm hook}(y h)
                                      -\rho xK_{\rm std}h\bigr).
\]

For the second identity differentiate $K_{\rm std}h$ in $x$.
Its linear factors are $y$ in the same-facet term and $D_ly$ in the
adjacent terms. Subtracting the first identity multiplied by $\rho$
leaves $D_l-C_l=w_lw_l^T$ and $w_l\cdot y=y_2-y_l$, proving
the displayed residual. Polynomial growth justifies all differentiations
and integrations by parts by Gaussian domination.

Setting $h=1$ in the first identity and using
$K_{\rm std}1=a_\rho/\rho$ recovers
$\mathsf K_{\rm hook}y=a_\rho x$, the rotation ground-state
identity, without a maximizing assumption. The second one gives
$\nabla a_\rho=(\rho^2/\sigma^2)\sum_lw_l\mathcal R_l1$.
Consequently

\[
 \nabla\bigl((a_\rho-K_{\rm std})h\bigr)
  =(a_\rho I-\rho\mathsf K_{\rm hook})\nabla h
      -\frac{\rho}{\sigma^2}\sum_lw_l
                     \bigl(\mathcal R_lh-\rho h\mathcal R_l1\bigr).
\]

Thus the kernels do not satisfy the plain OU derivative intertwining:
the explicit adjacent-facet residual must be controlled. Its displayed
identity does not assert a favorable sign for general balanced $h$,
and does not prove (S51). The exact checker controls the tangent-inverse
identity and the boundary-pair cancellation for $3\le k\le8$; the
unrestricted analytic identities are the written argument above.

## An exterior Gaussian identity settles the next complete hook

The junction cancellation in (S71) extends to every exterior degree.
It gives a strict supersolution for a further whole symmetry block,
including all its multiplicities.

### The complete exterior-hook kernels and their creation identity

For $0\le p\le k-2$ let $\mathsf K_p$ act on
$G_0$-equivariant $\Lambda^pH$-valued fields on $F_{12}$, with
$\mathsf K_0=K_{\rm std}$ and $\mathsf K_1=\mathsf K_{\rm hook}$.
The block is $\Lambda^{p+1}E\simeq V_{(k-p-1,1^{p+1})}$, and
its complete multiplicity space on $W$ is $L^2(W;\Lambda^pH)$.
The representation identification follows from the same vertical
Pieri input as (S4): the oriented $r$-subset module
$\Lambda^r\mathbb R^k$ has the two adjacent hook constituents,
whereas $\Lambda^r(\mathbf1\oplus E)=
\Lambda^rE\oplus\Lambda^{r-1}E$. Induction in $r$ gives the claim.
The $(12)$-minus space is $n\wedge\Lambda^pH$.

With $M_0,M_l$ as in (S71), its full-facet kernel is

\[
 \mathsf K_p(x,y)=2M_0(x,y)I+
                       \sum_{l=3}^k M_l(x,y)\Lambda^pC_l.
\]

Indeed the adjacent coefficient is

\[
 (2l)(n\wedge\omega)+(1l)(n\wedge\omega)
                         =n\wedge\Lambda^pC_l\omega.
\]

To verify it for all $p$, decompose $\omega=w_l\wedge\beta+\gamma$
with both $\beta,\gamma$ in exterior powers of $w_l^\perp\cap H$.
Both transpositions fix this orthogonal complement. On the plane
spanned by $n,w_l$, their tangent diagonal entries are $1/2,-1/2$,
their off-diagonal entries have opposite signs, and each determinant
is $-1$. Their sum therefore acts by $1$ on $n\wedge\gamma$ and
by $-2$ on $n\wedge w_l\wedge\beta$, exactly as the displayed
right side.

For a $G_0$-equivariant $C^1$ field
$\omega:H\to\Lambda^{p-1}H$ with polynomially bounded coefficients
and first derivatives, let $d\omega$ denote the exterior derivative.
Then, for $1\le p\le k-2$,

\[
 \boxed{\displaystyle
 \mathsf K_p(y\wedge\omega(y))(x)
    =\sigma^2\mathsf K_p(d\omega)(x)
                            +\rho x\wedge\mathsf K_{p-1}\omega(x).}
 \tag{S72}
\]

Apply Gaussian integration by parts to the exterior derivative.
The interior terms have the displayed form because $C_lD_l=I$:
applying $\Lambda^pC_l$ to $D_lx\wedge\omega$ gives
$x\wedge\Lambda^{p-1}C_l\omega$. To check every boundary term,
on the face $l$ the factor $2I+\Lambda^pC_l$ annihilates
$w_l\wedge\omega$, so the self and $l$th terms cancel.
For $r\ne l$ put $s=(lr)$ and $v=w_l-w_r$. Direct rank-one algebra
gives

\[
 C_rw_l=v,\qquad C_lw_r=-v,\qquad
 C_r-C_ls=v(w_l+w_r)^T.
\]

After the face substitution used in (S71), equivariance of $\omega$
makes the paired contribution proportional to
$v\wedge\bigl(\Lambda^{p-1}C_r-
\Lambda^{p-1}(C_ls)\bigr)\omega=0$.
The last equality holds because the two linear maps differ by a
map with range spanned by $v$; expanding their exterior powers
leaves a factor $v$ in every difference term. Gaussian decay removes
the boundary at infinity. This proves (S72) with the stated regularity.

The folded kernels preserve a useful fixed coordinate cone in every
exterior degree. Define the orthonormal Helmert vectors

\[
 u_j=\frac{e_1+\cdots+e_j-je_{j+1}}{\sqrt{j(j+1)}},
                       \qquad1\le j\le k-1.
\]

Here $u_1=n$ and $u_2,\ldots,u_{k-1}$ form a basis of $H$.
For $j\ge2$, the adjacent transposition $(j,j+1)$ mixes only
$u_{j-1},u_j$, with positive
off-diagonal entries. In the increasing exterior basis its off-diagonal
entries are still positive: the mixed indices are consecutive, so
reordering introduces no extra sign. A wedge containing both mixed
indices instead has eigenvalue $-1$. Thus every $I+R(s_j)$ is
entrywise nonnegative in this exterior basis, and $(12)$ is diagonal.
The interval-symmetrizer argument of (S69), applied in this same
basis and then compressed to wedges containing $u_1$, proves that
the folded $\mathsf K_p$ is entrywise nonnegative in the increasing
basis of $\Lambda^pH$. This uses the retained total-positivity input;
no general spectral ordering between different hooks is assumed.

### A strict supersolution for the third exterior hook

**Theorem S73 (third exterior hook).**

For every $k\ge4$, $d\ge k-1$ and $0<\rho<1$, the entire
$V_{(k-3,1,1,1)}$ isotypic component of the cylindrical simplex interface form has
a strictly positive $L^2$ gap and no kernel:

\[
 \mathfrak q_\rho(f)\ge\kappa^{(3)}_{k,\rho}\|f\|^2,
                    \qquad\kappa^{(3)}_{k,\rho}>0.
 \tag{S73}
\]

The constant is qualitative here and independent of $d$. The assertion includes all fields
in the multiplicity space, not only polynomial or radial speeds.
For $k=4$ this is another proof of the sign gap; for $k\ge5$ it
settles an additional hook beyond (S56).

Put $p_0=(e_1+e_2)/2-k^{-1}\mathbf1\in H$, and define
$W_2(x)=p_0\wedge x$. Write
$R_l(x)=\int_{F_{12}}M_l(x,y)d\nu_{12}(y)>0$ and
$v(x)=\sum_{l=3}^kR_l(x)w_l$. Since
$w_l\cdot p_0=1/2$, (S49) gives

\[
 \mathsf K_1p_0=(a_\rho/\rho)p_0-v.
\]

Use (S72) at $p=2$ with the constant one-form $\omega=p_0$.
Its exterior derivative is zero. Reversing the wedge order yields
the exact supersolution residual

\[
 \mathsf K_2W_2=a_\rho W_2-\rho v\wedge x.
\]

Both $W_2$ and $v\wedge x$ have strictly positive coordinates in
the increasing Helmert exterior basis on $W$. Here are the required
inequalities in full. Set $\alpha_j=1/\sqrt{j(j+1)}$ and

\[
 X_j=\frac{\langle x,u_j\rangle}{\alpha_j}
       =\sum_{i=1}^jx_i-jx_{j+1},\qquad
 V_j=\frac{\langle v,u_j\rangle}{\alpha_j}
       =(j+1)R_{j+1}+\sum_{l=j+2}^kR_l,
                       \quad2\le j\le k-1.
\]

On $W$, $0<X_2<X_3<\cdots<X_{k-1}$. Also
$R_3>R_4>\cdots>R_k>0$: to compare $R_j,R_{j+1}$ substitute
$y\mapsto(j,j+1)y$ in the latter integral. The resulting exponent
difference is
$\tau(x_j-x_{j+1})(y_2-y_{j+1})>0$ on the facet interior.
Consequently

\[
 V_{j+1}-V_j=(j+1)(R_{j+2}-R_{j+1})<0,
                   \qquad V_j>0.
\]

Since $\langle p_0,u_j\rangle=\alpha_j$, the $i,j$ coordinates
for $2\le i<j\le k-1$ are

\[
 (W_2)_{ij}=\alpha_i\alpha_j(X_j-X_i)>0,\qquad
 (v\wedge x)_{ij}
      =\alpha_i\alpha_j(V_iX_j-V_jX_i)>0.
\]

The last strict sign follows from $V_i>V_j>0$ and $X_j>X_i>0$.
Thus, with the folded operator,
$T_2=a_\rho^{-1/2}\mathsf K_2a_\rho^{-1/2}$ has the positive
$L^2$ supersolution $\psi=\sqrt{a_\rho}W_2$ satisfying
$T_2\psi<\psi$ in every coordinate almost everywhere.
Its kernel is entrywise nonnegative in these orthonormal coordinates,
and it is compact and self-adjoint by the full interface bounds.
For any eigenvector $z\ne0$, integration of
$|T_2z|\le T_2|z|$ against $\psi$ gives

\[
 |\lambda|\langle\psi,|z|\rangle
       \le\langle T_2\psi,|z|\rangle
           <\langle\psi,|z|\rangle.
\]

Hence every eigenvalue has modulus less than $1$, and compactness
gives $\theta_{k,\rho}=\|T_2\|<1$. Taking
$\kappa^{(3)}_{k,\rho}=\mu_\rho(1-\theta_{k,\rho})>0$ proves
the active-space case for the full Hilbert-space block. No mass
constraint occurs in this representation.

For $d>k-1$, the product Mehler kernel and Gaussian surface measure
factor exactly as in (S38). Expanding in an orthonormal Hermite basis
in the unused coordinates multiplies the interaction at degree $m$
by $\rho^m$, leaving $a_\rho$ unchanged. Its normalized norm is
$\rho^m\theta_{k,\rho}\le\theta_{k,\rho}$, so every degree has
the same lower gap bound. Parseval proves (S73) for all $L^2$ speeds
in every ambient dimension, with the stated dimension-independent
constant.

The third-hook result alone leaves the standard $(4,1)$ and
$(2,2,1)$ blocks in five cells. The next section settles the latter;
(S73) does not settle every higher exterior degree.

## Five cells: an exact positive polynomial certificate for the mixed block

**Theorem S74 (five-cell mixed block).**

For $k=5$, every $d\ge4$ and every $0<\rho<1$, the entire
$V_{(2,2,1)}$ isotypic component has a positive $L^2$ gap and no
kernel:

\[
 \mathfrak q_\rho(f)\ge\kappa^{221}_{5,\rho}\|f\|^2,
               \qquad\kappa^{221}_{5,\rho}>0.
 \tag{S74}
\]

The constant is independent of $d$ and qualitative here. The proof
uses an exhaustive rational polynomial expansion, with a written
covering argument for every strictly totally positive matrix. It is
not an extrapolation from selected chamber points or correlations.
The analytic operator proof and the finite calculation have only
internal review; there is no proof-assistant formalization or
independent replay.

### The linear positive field and exact residual

In a rational seminormal basis of the $(12)$-minus space, put

\[
 L=\begin{pmatrix}
 0&0&0&5/4&-5/4\\
 0&0&5/4&-5/12&-5/6\\
 0&0&5/8&-5/8&0
 \end{pmatrix},\qquad W(x)=Lx.
\]

This basis is a positive diagonal rescaling of the orthonormal Young
basis, with metric $J=\operatorname{diag}(1,9/8,3/2)$.
The adjacent transports $T_l=2P_-R((2l))P_-$ are

\[
 T_3=\begin{pmatrix}-2&0&0\\0&1&0\\0&0&1\end{pmatrix},\quad
 T_4=\begin{pmatrix}2/3&-1&0\\-8/9&-5/3&0\\0&0&1\end{pmatrix},\quad
 T_5=\begin{pmatrix}2/3&1/2&1\\4/9&1/3&-4/3\\2/3&-1&-1\end{pmatrix}.
\]

The checker reconstructs these data from Young's rational generators
for shape $(2,2,1)$, verifies all Coxeter relations and the
$G_0=S_{\{3,4,5\}}$ equivariance $L\pi=R(\pi)L$.
To specify the basis without a phase convention, the three retained
tableaux have the box positions of labels $1,\ldots,5$ respectively

\[
 ((0,0),(1,0),(2,0),(0,1),(1,1)),\quad
 ((0,0),(1,0),(0,1),(2,0),(1,1)),\quad
 ((0,0),(1,0),(0,1),(1,1),(2,0)).
\]

The rational generators have diagonal $1/r$ at axial distance $r$;
the off-diagonal entry towards the tableau with longer row-reading
inversion count is $1$, and the reverse entry is $1-r^{-2}$.
These explicit finite matrices, rather than floating point square
roots, are the inputs to the coefficient calculation.

On $W$, all components of $W(x)$ are positive: they are positive
multiples of $x_4-x_5$, $3(x_3-x_4)+2(x_4-x_5)$ and
$x_3-x_4$. Also $Lp_0=0$. The full-facet kernel is
$\mathsf K_{221}=2M_0I+\sum_{l=3}^5M_lT_l$, and its folded
kernel is entrywise nonnegative by (S67)--(S69). Its self-adjoint
metric is $J$, as verified by $JT_l=T_l^TJ$ and the corresponding
unitarity of $G_0$.

Let $m_l(x)=\int_{F_{12}}M_l(x,y)y\,d\nu_{12}(y)$, including
$l=0$. The rotation identity (S71) says
$2m_0+\sum_lC_lm_l=a_\rho x$. Applying $L$ and subtracting the
mixed-block interaction gives the exact identity

\[
       a_\rho W-\mathsf K_{221}W
                =\sum_{l=3}^5(LC_l-T_lL)m_l.
\]

It remains to prove that this residual is strictly positive in every
coordinate. That is the purpose of the following finite formal
calculation.

### Nine immanants positive on every strictly TP four-by-four matrix

Let $r_m$ be the chamber rays of (S53), now with $k=5$, and
$y=\sum_{m=2}^4\eta_mr_m$, where
$\eta_m=(y_m-y_{m+1})/5>0$. For $\pi\in G_0$ and $3\le l\le5$
put $g=(2l)\pi\in S_{\{2,3,4,5\}}$. Define a rational coefficient
vector indexed by this group by

\[
 c_{g,m}=(LC_l-T_lL)\pi r_m,
\]

and set it to zero for permutations fixing $2$. The displayed
parametrization of the other permutations is unique. For a $4\times4$
matrix $A$ indexed by $2,3,4,5$, define the three-component immanant

\[
 \Phi_m(A)=\sum_g c_{g,m}\prod_{i=2}^5A_{i,g^{-1}(i)}.
 \qquad (\Phi_m(A))_a>0
       \quad\text{for every strictly TP }A,\ a=1,2,3.
 \tag{S75}
\]

Here is a covering reduction that makes a finite coefficient check
sufficient for the stated unrestricted matrix assertion. Write
the four matrix indices as $1,2,3,4$ in this factorization. Set
$x_i(t)=I+tE_{i+1,i}$ and $z_i(t)=I+tE_{i,i+1}$, and substitute

\[
 A=x_1(t_1)x_2(t_2)x_1(t_3)x_3(t_4)x_2(t_5)x_1(t_6)
       \operatorname{diag}(d_1,d_2,d_3,d_4)
       z_1(s_1)z_2(s_2)z_3(s_3)z_1(s_4)z_2(s_5)z_1(s_6),
\]

with all sixteen parameters positive. Every strictly TP $4\times4$
matrix has this form. For a direct verification, ordinary elimination
gives $A=BDU$ with unit lower and upper factors and positive diagonal
pivots. For each $r$, the flag minors of $B$ with columns $1,\ldots,r$
are the corresponding minors of $A$ divided by $d_1\cdots d_r$,
so they are positive. Write

\[
 B=\begin{pmatrix}1&0&0&0\\a&1&0&0\\b&c&1&0\\d&e&f&1\end{pmatrix}.
\]

Its lower Jacobi parameters, in the order used above, are

\[
 t_1=\frac{a(cf-e)-bf+d}{cf-e},\quad
 t_2=\frac{cf-e}{f},\quad
 t_3=\frac{f(be-cd)}{e(cf-e)},\quad
 t_4=f,\quad t_5=e/f,\quad t_6=d/e.
\]

Every numerator and denominator is positive by those flag minors:
$d,e,f,cf-e,be-cd,a(cf-e)-bf+d$ are respectively the flag
minors with row/column sets
$(4;1)$, $(1,4;1,2)$, $(1,2,4;1,2,3)$,
$(1,3,4;1,2,3)$, $(3,4;1,2)$ and $(2,3,4;1,2,3)$.
Multiplying the six elementary factors gives $B$ with these parameters.
Applying the same argument to $A^T$ supplies the positive parameters
for $U$ in the reversed upper word. Thus no strictly TP input is
excluded by the substitution.

The exact checker constructs the sixteen-variable matrix and all
twenty-four permutation monomials, derives the $c_{g,m}$ from the
displayed representation data, and collects every coefficient using
rational arithmetic. The numbers of nonzero monomials are

| Component | $m=2$ | $m=3$ | $m=4$ |
| --- | ---: | ---: | ---: |
| $1$ | 164 | 158 | 187 |
| $2$ | 290 | 224 | 217 |
| $3$ | 126 | 78 | 61 |

All **1505** coefficients are strictly positive. The receipt records
the count, minimum coefficient, maximum total degree and a canonical
coefficient hash for each polynomial. The whole polynomial is
reconstructed from the retained code and inputs, with no interpolation
or numerical tests. Each polynomial is therefore strictly positive
for all positive parameter values. The covering argument proves (S75).

Total positivity is necessary for this method. For the positive-entry
matrix with rows $(3,5,5,3)$, $(2,6,6,5)$, $(4,4,4,1)$ and
$(6,3,6,6)$, the minor with rows $1,3$ and columns $1,2$ is $-8$,
and the second component of $\Phi_2$ is $-22825$. The checker verifies
this rejection control. It is not a Gaussian simplex counterexample.

### Strict supersolution and the five-cell reduction

For $x,y\in W$, use $A=[e^{\tau x_i y_j}]_{i,j=2}^5$, which is
strictly TP by (S47). Splitting the residual integral into its six
remaining-coordinate chambers gives

\[
 a_\rho(x)W(x)-(\mathsf K_{221}W)(x)
  =\int_W\Gamma_\rho(x,y)e^{\tau x_2y_2}
                    \sum_{m=2}^4\eta_m\Phi_m(A)\,d\nu_{12}(y)>0
\]

in every coordinate. The positive field $\psi=\sqrt{a_\rho}W$
is consequently a strict supersolution for the normalized folded
operator. Rescale by $J^{1/2}$ to make the basis orthonormal; this
preserves all coordinate inequalities and kernel nonnegativity.
The eigenvector integration and compactness argument used for (S73)
give norm $\theta_{221,\rho}<1$ and active gap
$\mu_\rho(1-\theta_{221,\rho})>0$ on the entire multiplicity space.
The same Hermite factorization as in (S73) proves (S74) in every
ambient dimension, with this unchanged gap. The positive field is a
device for bounding the whole operator, not a restriction to linear
speeds.

The nontrivial $S_5$ types are $(4,1)$, $(3,2)$, $(3,1,1)$,
$(2,2,1)$, $(2,1,1,1)$ and $(1^5)$. The moment block is covered by
(S60), the rotation hook by (S56), the mixed block by (S74), the
third exterior hook by (S73), and the sign by (S62). Hence for each
fixed $0<\rho<1$,

\[
 \boxed{\text{Full five-cell S37-type gap modulo rotations}
                     \quad\Longleftrightarrow\quad
                       \text{the five-cell standard inequality (S51).}}
 \tag{S76}
\]

For sufficiency, (S51) has a positive gap by the compact normalized
operator and the positive lower bound on $a_\rho$. Taking the minimum
over the six symmetry types gives the active gap, and (S42) lifts it
uniformly in dimension. Necessity follows by restricting to the active
standard block. This equivalence does not prove the remaining standard
inequality, which retains the precise two-part test (S70).

## The standard block at small correlation and a limit of the even comparison

Two further facts distinguish the higher-cell standard problem from a
direct repetition of (S32). Its initial interaction has a simple exact
rank-one description, but the elementary antipodal pointwise comparison
used in the three-cell proof fails in higher cell counts.

### The exact initial standard spectral minimum

Put $v_k=\nu_{12}(F_{12})$, $\mu=\nu_{12}/v_k$ on the full
facet, $p=p_0/|p_0|$, and $s_k^2=\operatorname{Var}_\mu(p\cdot x)$.
For $G_0$-invariant $h$ with $\int h\,d\nu_{12}=0$, define
$q_{{\rm std},\rho}(h)=\rho^{-1}
(\int a_\rho h^2d\nu_{12}-\iint K_{{\rm std},\rho}hh\,d\nu_{12}^2)$.
Then for every fixed $k\ge3$,

\[
 \boxed{\displaystyle
 \lim_{\rho\downarrow0}\ inf_{\substack{\|h\|_{L^2(\nu_{12})}=1\\
           h\ G_0\text{-invariant},\ \int h\,d\nu_{12}=0}}
 q_{{\rm std},\rho}(h)
       =\ell\left(1-\frac{s_k^2}{2}\right)\ge\frac\ell2.}
 \tag{S77}
\]

In particular this block has gap at least $\ell/3$ on some
fixed-$k$ positive-correlation interval, without using the accepted
full first-moment stability theorem. No interval endpoint is evaluated
here and no conclusion for every correlation follows from this limit.

First $\ell=kv_k$: Gaussian divergence on $T_i$ gives
$b_i^*=\sum_{j\ne i}v_kn_{ij}=(kv_k/\sqrt2)(e_i-k^{-1}\mathbf1)$.
The difference of two such moments has length $kv_k$.
The row identity (S49) gives $a_\rho/\rho\to kv_k$ pointwise.
On the mean-zero scalar space, the compressed interaction
$K_{{\rm std},\rho}/\rho$ converges in Hilbert--Schmidt norm to the
compressed derivative of its kernel at zero. This derivative is

\[
        2x\cdot y+\sum_{l=3}^k x\cdot(2l)y
                  =x\cdot\left(kI-\sum_lw_lw_l^T\right)y.
\]

Indeed the constant kernel is $k$ and is annihilated by the mean-zero
compression. For $|\rho|\le1/4$, two derivatives of the kernels
have a common integrable Gaussian-times-polynomial square bound,
so Taylor's formula proves the claimed Hilbert--Schmidt convergence
with $O_k(\rho)$ error after division by $\rho$.
Symmetry forces $m_h=\int xh\,d\nu_{12}$ into $\mathbb Rp_0$.
Since $\sum_lw_l=kp_0$ and $w_l\cdot p_0=1/2$, the limiting
interaction is exactly $(k/2)|m_h|^2$. Its only nonzero eigenvalue
on this scalar space is $(k/2)v_ks_k^2=\ell s_k^2/2$, attained
by the centered function $p\cdot x-\mathbb E_\mu[p\cdot x]$.

The classical convex-Gaussian Poincare bound gives $s_k^2\le1$.
The precise inspected input is Brandolini--Chiacchio--Henrot--Trombetti,
Theorem 1.1 and the mean-zero Rayleigh characterization, pp. 1--2 of
the [retained arXiv v2](sources/brandolini-chiacchio-henrot-trombetti-poincare-12096469v2.pdf).
It concerns convex $C^2$ domains with a uniform interior sphere
condition; a polyhedral facet is not directly within those hypotheses.
For the application here, approximate each bounded convex set
$F_{12}\cap B_R\subset H$ by bounded smooth strictly convex domains.
Apply that theorem to the linear function centered by its conditional
Gaussian mean. Its gradient has norm $1$, so its conditional variance
is at most $1$. Gaussian moments of orders zero through two converge
under the bounded convex approximation. Then let $R\to\infty$;
Gaussian moment domination proves the same variance bound on $F_{12}$.
The sphere radius may depend on the approximating domain; the
Poincare constant $1$ does not. For $k=3$, the facet is a halfline
and $s_3^2=1-2/\pi\le1$ directly. We use the mean-zero Rayleigh
statement, with the conditional mean normalized by the domain's
Gaussian measure.

Finally (B7) gives $a_\rho/\rho\ge\ell/\sigma$. Together with
the interaction norm limit this proves the uniform lower limit in
(S77). Testing its centered linear eigenfunction gives the matching
upper limit by dominated convergence, using the upper bound in (B7).
This argument does not assert uniform convergence of the multiplier
on the unbounded facet. That convergence is unnecessary and generally
fails. The covariance bound is an imported analytic input, not a
consequence of the finite checker.

### A pointwise antipodal OU comparison fails from seven cells onward

Let the candidate antipodal folded kernel be

\[
 \mathcal E_\rho(x,y)=\Gamma_\rho(x,y)
       \sum_{\pi\in G_0}
          \left(e^{\tau x\cdot\pi y}+e^{-\tau x\cdot\pi y}\right).
\]

It is the even OU kernel on the supporting space $H$, folded by
$G_0$ and restricted to the facet, up to a fixed positive normalization.
For every fixed $k\ge7$ and $0<\rho<1$,

\[
 \boxed{\displaystyle
  \inf_{x,y\in W}
       \frac{\mathcal K_{{\rm std},\rho}(x,y)}{\mathcal E_\rho(x,y)}=0.}
 \tag{S78}
\]

Thus no positive constant depending only on $k,\rho$ supplies this
particular pointwise lower comparison. This is not a counterexample
to (S51), nor does it rule out an integrated or different reflected
comparison.

For the proof begin with chamber rays $x=r_3$, $y=r_{k-1}$.
Rearrangement gives $x\cdot gy\le x\cdot y=3k$ for every
$g\in S_{\{2,\ldots,k\}}$. The remaining-label transposition
$\pi_*=(3k)$ instead gives
$x\cdot\pi_*y=k(3-k)$. Hence
$-x\cdot\pi_*y-x\cdot y=k(k-6)>0$.
This strict inequality persists at interior chamber points. For an
explicit choice put $R=\sum_{m=2}^{k-1}r_m$,
$\epsilon=k^{-4}$, $x_\epsilon=r_3+\epsilon R$ and
$y_\epsilon=r_{k-1}+\epsilon R$. Both lie in $W$. The perturbation
of the strict difference has magnitude at most
$3/\sqrt k+1/(2k^3)<2$, using
$|R|\le k^{5/2}/2$ and $|r_3|+|r_{k-1}|<3k$; the original
difference is at least $7$.

Set $A=x_\epsilon\cdot y_\epsilon$ and
$B=-x_\epsilon\cdot\pi_*y_\epsilon>A$. The folded standard
kernel has $k(k-2)!$ positive exponential terms counting its
coefficients, all with exponent at most $\tau A$. The antipodal
kernel contains the term with exponent $\tau B$. At scaled points
$t x_\epsilon,t y_\epsilon$ their common $\Gamma_\rho$ cancels,
leaving

\[
 \frac{\mathcal K_{{\rm std},\rho}(t x_\epsilon,t y_\epsilon)}
      {\mathcal E_\rho(t x_\epsilon,t y_\epsilon)}
   \le k(k-2)!\exp\!\left(-\tau t^2(B-A)\right)\longrightarrow0.
\]

This proves (S78) at every fixed positive correlation. The checker
verifies the ray algebra, the exact interior perturbation and exponent
separation for $7\le k\le20$, with a $k=6$ equality control. The
unbounded-scale conclusion follows from the displayed exponential
bound, not from sampled values of the kernel.

## Convex Gaussian Poincare does not give the missing nonlocal mass estimate

The local covariance argument in (S77) cannot simply be extended to a
nonlocal comparison for every convex Gaussian domain. In fact the
following candidate inequality is false:

\[
 \iint_{\Omega^2}M_\rho(x,y)h(x)h(y)d\gamma_d(x)d\gamma_d(y)
       \le\rho\int_\Omega R_\rho(x)h(x)^2d\gamma_d(x),\qquad
 R_\rho(x)=\int_\Omega M_\rho(x,y)d\gamma_d(y),\quad
 \int_\Omega h\,d\gamma_d=0.
\]

It would have the right ordinary-mass constraint and row multiplier
for a one-step Poincare comparison. However, for **every** fixed
$0<\rho<1$, it fails on a positive orthant in some finite dimension,
even for a permutation-invariant linear $h$.
At $\rho=1/2$ there is the explicit counterexample

\[
 \boxed{\begin{gathered}
 \Omega=(0,\infty)^{19},\qquad
 h(x)=\sum_{i=1}^{19}x_i-19\sqrt{\frac2\pi},\\
 \frac12\int_\Omega R_{1/2}h^2d\gamma_{19}
     -\iint_{\Omega^2}M_{1/2}h(x)h(y)d\gamma_{19}^2
       =\left(\frac13\right)^{19}
                 \frac{19(62-36\sqrt3)}{64\pi}
                    <-\frac1{32}\left(\frac13\right)^{19}<0.
 \end{gathered}}
 \tag{S79}
\]

This does not refute (S51): the simplex standard kernel includes the
specific adjacent-facet terms, and its domain is not this orthant.
It rules out deriving the needed nonlocal inequality from convexity,
Gaussianity and the ordinary mean constraint alone.

Here is a self-contained calculation. Let $(X,Y)$ be a standard
normal pair with correlation $\rho$, and write
$A=\mathbb P(X>0,Y>0)=1/4+\arcsin(\rho)/(2\pi)$.
Put $\mu_0=\mathbb E[X\mid X>0]=\sqrt{2/\pi}$ and
$\alpha=\mathbb E[X\mid X>0,Y>0]$. Gaussian integration by parts
gives

\[
 \alpha=\frac{(1+\rho)\phi(0)}{2A},\qquad
 D=\mathbb E[X^2\mid X,Y>0]
       =1+\frac{\rho\sqrt{1-\rho^2}}{2\pi A},\qquad
 C=\mathbb E[XY\mid X,Y>0]
       =\rho+\frac{\sqrt{1-\rho^2}}{2\pi A}.
\]

For example, write $c=\rho/\sqrt{1-\rho^2}$ and integrate
$x\phi(x)\Phi(cx)$ and $x^2\phi(x)\Phi(cx)$ on the halfline.
For the cross moment use
$\mathbb E[Y1_{\{Y>0\}}\mid X=x]
=\rho x\Phi(cx)+\sqrt{1-\rho^2}\phi(cx)$.
The quadrant probability follows from the Gaussian polar angle.
These formulas use the Gaussian-relative Mehler kernel and no surface
measure normalization.

Put $h_d(x)=\sum_i x_i-d\mu_0$. Its ordinary Gaussian integral
on $(0,\infty)^d$ is zero by independence, and it belongs to $L^2$.
Under the conditional joint law of all $d$ pairs in the positive
quadrant, its mean is $d(\alpha-\mu_0)$, its variance is
$d(D-\alpha^2)$ and its cross covariance is $d(C-\alpha^2)$.
Consequently the exact normalized difference is

\[
 \frac{\rho\int R_\rho h_d^2d\gamma_d
                  -\iint M_\rho h_dh_d\,d\gamma_d^2}{A^d}
   =d\bigl[\rho D-C+(1-\rho)\alpha^2\bigr]
                     -(1-\rho)d^2(\alpha-\mu_0)^2.
\]

For $0<\rho<1$, strict convexity of $\arcsin$ between its endpoints
gives $2\arcsin(\rho)/\pi<\rho$, hence $\alpha>\mu_0$.
The negative term is therefore quadratic in $d$ with a strictly
negative coefficient, while the other term is linear. This proves the
all-correlation counterexample family in finite, sufficiently large
dimensions. It identifies the mass-weight mismatch concretely: the
ordinary centering is $d\mu_0$, whereas the pair-conditioned measure
has mean $d\alpha$.

For the explicit instance in (S79), $A=1/3$,
$\alpha=(9/8)\mu_0$, $D=1+3\sqrt3/(8\pi)$ and
$C=1/2+3\sqrt3/(4\pi)$. Substitution simplifies the difference in
dimension $d$ to

\[
       \frac{d(81-d-36\sqrt3)}{64\pi}.
\]

At $d=19$, its sign is negative since $62^2<36^2\cdot3$.
For the quantitative margin use $\sqrt3>433/250$ and
$\pi<22/7$, the latter proved in (S36). Since the numerator is
negative, directed arithmetic gives

\[
       \frac{19(62-36\sqrt3)}{64\pi}
           <-\frac{133}{4000}<-\frac1{32}.
\]

This proves the bound in (S79). In this particular linear family the
first negative integer dimension at $\rho=1/2$ is $19$:
$d=18$ is positive since $63^2>36^2\cdot3$. In dimension one the
normalized difference is $(20-9\sqrt3)/(16\pi)>0$.
This is a statement about the displayed family, not a classification
of all counterexamples in lower dimensions. The ordinary-centering
estimate does not tensorize in the way a local Poincare inequality does.

Indeed the local Gaussian Poincare gap on an orthant is $2$, by
even extension in every coordinate and the Gaussian Hermite spectrum.
Even that sharper local fact is consistent with this nonlocal
counterexample. The checker certifies the exact quadratic-in-dimension
moment identity, the rational negative margin and both lower-dimensional
sign controls; the pair moment integrals, orthant realization and
all-correlation family use the written proof.

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
or a different quantity. The
[perimeter-angle calculation](THREE_CELL_DEFICIT_AND_PERIMETER.md)
also refutes the universal nondecreasing law for
$\Delta_\rho/\arccos\rho$ by a balanced three-slab partition.

The remaining targets are distinct:

* Prove the active-space higher-cell gap in (S4). The standard block
  has the positive kernel (S49)--(S50)
  and the unresolved weighted Poincare inequality (S51).
  Its exact two-part spectral criterion is (S70); a bound on the
  second normalized eigenvalue alone does not settle ordinary balance.
  Its small-correlation spectral limit is (S77). For $k\ge7$, (S78)
  rules out directly extending the pointwise antipodal comparison (S32).
  Equation (S79) also excludes a generic convex-domain ordinary-mass
  nonlocal Mehler inequality; the special adjacent-facet structure
  must enter any use of that route.
  The complete rotation hook has the all-correlation gap (S56).
  The $(k-2,2)$ block has the strict gap (S60), and every sign block is
  settled by (S61)--(S66). The third exterior hook is settled by (S73).
  In four and five cells only (S51) remains, by (S76) and its four-cell
  analogue. For $k\ge6$ the other blocks also need estimates.
  Equation (S42) then handles every ambient dimension.
  Equation (S23) and the full three-cell interface problem are settled
  by (S29)--(S38).
* Classify stationary stable balanced partitions. Simplex spectral
  stability does not classify critical points away from the simplex.
* Exclude a nonsimplex maximizer tied with the simplex. The continuation
  dichotomy retains this alternative even if the entire simplex Hessian
  is strictly positive modulo rotations.
* Prove exact optimality near $\rho=1$ using a quantitative perimeter
  bound together with estimates relative to the deformation of a moving
  competitor. A perimeter limit for each fixed partition is insufficient.
* Seek stronger moment or perimeter inequalities that control that
  relative error, and test a modified monotonicity statement that
  accounts for the simplex junction term.

Run the exact checker and its optional exploratory quadrature with

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python research/plurality-stability/check_simplex_interface_symmetry.py --quadrature --output research/plurality-stability/simplex_interface_symmetry_verification.json
```

The exact part requires only Python's standard library. It checks all 30
permutations in the $S_3/S_4$ character examples, (S21), the elementary
polynomial bounds (S32)/(S36), and the interval inequalities used in
(S27)--(S28), including rejection controls. It also compares the geometric
coset sum (S44) with the determinant (S46) as exact formal exponential
sums for $3\le k\le8$, including omitted-term and wrong-sign controls,
and verifies the standard-sector permanent identity (S49)--(S50).
It checks the adjacent bivector identity and chamber Gram matrix for
$3\le k\le7$, and exact positive folded hook matrices at one stated
pair of chamber points for each such $k$, with $\tau=\log2$ so every
exponential is an exact power of two. It rejects the wrong reflection
coefficient and a point outside the ordered cone. These are bounded
controls; it also checks the symmetric transport, two-row Gram and
ordered-residual exponent identities, selected positive two-row matrices,
and one exact four-cell sign domination. The universal pairing,
supersolution and spectral arguments are written proofs.
For (S63)--(S66) it checks selected nonsymmetric exponential-matrix
minor inequalities and exact Gaussian Vandermonde moments for one
through six untied coordinates, with controls for the wrong correlation
degree, insufficient positivity hypotheses and the $k=3,4$ exceptional
scopes. It does not verify the imported Hermitian transition theorem.
For (S67)--(S69) it checks the free group-algebra identity through $k=7$,
exact rational Young-generator Coxeter relations and selected positive
immanant matrices for every partition through $k=5$ and four partitions
at $k=6$. The unrestricted matrix-positivity theorem uses the two retained
primary statements and the written compression argument.
It also verifies the exact two-atom mass-weight counterexample following
(S70). The infinite-dimensional equivalence (S70) is a written proof,
and neither of its two tests is certified for the general simplex.
For (S71) it checks the tangent-compression inverse and all boundary-pair
permutation identities through eight cells, not the analytic integrations
by parts or a sign bound for their residual.
For (S72)--(S73) it checks adjacent third-exterior transport and the
rank-one boundary-pair identity through eight cells, selected positive
folded exterior matrices, and exact increasing/decreasing coordinate
controls. The Gaussian creation identity and the all-correlation
supersolution proof remain written analytic arguments.
For (S75) the checker exhaustively expands nine fixed rational
polynomials in sixteen formal variables and checks all 1505 nonzero
coefficients. The unrestricted matrix conclusion also requires the
written positive-factorization covering argument. The passage from
these polynomials to the full five-cell block gap and (S76) is analytic
and representation-theoretic written mathematics.
For (S77)--(S78) it checks the rank-one transport algebra, ray dot
products and exact interior perturbations in stated finite cell ranges.
The Gaussian covariance input, the Hilbert--Schmidt Taylor limit and
the unrestricted exponential-scale obstruction are written arguments.
For (S79) it checks the dimension polynomial and the explicit
19-dimensional normalized negative margin. The Gaussian pair moments,
orthant independence and all-correlation family are analytic
written arguments, not assertions inferred from a bounded dimension scan.
The analytic operator reduction, OU parity spectra, tensor decomposition,
dimension lifting, determinant positivity and all-$k$ representation
argument still require written proof review. The optional
quadrature uses NumPy 1.26.4 to check the displayed gap constants and
extra-degree-one block at four chosen correlations, with two grids on
$[0,10]$, and checks the separate perimeter expansion at three angles.
It has no certified tail or discretization errors and proves no
inequality. Its agreement is separate from the exact arithmetic.
