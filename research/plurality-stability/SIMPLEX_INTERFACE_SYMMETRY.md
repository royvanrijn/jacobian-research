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
The single-reflection argument (S85)--(S87) further settles every
nonstandard two-row block $(k-p,p)$, $2\le p\le\lfloor k/2\rfloor$,
at every positive correlation, with one gap constant valid for all $p$.
The antisymmetry hierarchy (S88)--(S90) identifies the next obstruction
and embeds every block of each fixed height into one scalar spectral
problem, with an explicit closed physical image. The highest sign
projection loses no equivariant field. These reductions supply no new
higher-row gap by themselves.
The exterior Gaussian identity (S72) and a linear positive
supersolution also settle the whole third exterior hook
$(k-3,1,1,1)$ for every $k\ge4$, by (S73).
Equation (S91) supplies a positive harmonic polynomial field in every
higher exterior degree and an exact one-dimensional formula for its
truncated image. The convex clipping argument (S92) settles the fourth
exterior hook $(k-4,1^4)$ for every $k\ge5$, with the explicit gap
$\ell\sqrt{1-\rho^2}$. The sufficient comparison remains open in the
other intermediate exterior degrees.
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
The Hermite floor (S80) and Schur estimate (S81) also give explicit
gaps in every sufficiently deep representation at each fixed correlation,
including the complete two-column family $(2^r)$, $r\ge4$, for
$0<\rho\le1/2$.
The convex covariance proof (S82) gives the explicit first-moment
interface gap $\ell/k$ without global first-moment optimality. It
also supplies a uniform initial full gap $\ell/(2k)$ in (S83), on
the explicit interval $0<\rho\le[6\sqrt3\,k(k-1)]^{-1}$ in
every ambient dimension.

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
with kernel $\mathcal R_k$. The independent convex-covariance proof
(S82) below gives the stronger explicit gap $\ell/k$ without the
global first-moment theorem; it does not use this paragraph.
Polarization at $r\in\mathcal R_k$ gives

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

## Every representation has an explicit Hermite floor

The sign cancellation (S12) is a special case of a representation-wide
degree obstruction. It gives useful bounds even when a positive
supersolution for the compressed immanant has not been found.

For a nontrivial partition $\lambda=(\lambda_1,\ldots,\lambda_r)$ of
$k$, put

\[
 n(\lambda)=\sum_{i=1}^r(i-1)\lambda_i,\qquad
 D_\lambda=n(\lambda)-(r-1).
\]

On its complete interface isotypic space, in every ambient dimension
$d\ge k-1$, the coefficients in (S11) satisfy

\[
 \boxed{E_m(f)=0\quad(m<D_\lambda),\qquad
 0\le B_\rho(f)\le(\rho/b)^{D_\lambda-1}B_b(f)
                  \quad(0<\rho\le b<1).}
 \tag{S80}
\]

The comparison includes $D_\lambda=0$, although only $D_\lambda\ge2$
will give the new gap estimates below. For the standard, rotation hook,
two-row moment block, mixed $(k-3,2,1)$ block and third exterior hook,
the respective floors are $0,1,1,2,3$. For the sign they are
$(k-1)(k-2)/2$, recovering (S12).

Here is the representation argument, including the polynomial input.
In characteristic zero the smallest polynomial degree of a Specht
module $V_\beta$ in independent variables permuted by $S_s$ is
$n(\beta)$. The Specht polynomial is the product of the Vandermondes
of the columns, of that degree. The lower bound follows because a
nonzero Young-symmetrizer image alternates within each column and is
therefore divisible by that product. For the precise classical facts,
see Debus--Kretschmer, *Symmetric Ideals and Invariant Hilbert Schemes*,
[arXiv:2404.15240v2](https://arxiv.org/pdf/2404.15240v2), section 2.1,
Theorem 2.1 and the paragraph after Lemma 3.7 on printed page 9.
This is a classical polynomial representation input, not a new
minimal-degree theorem here. The retained PDF was inspected, including
the column-divisibility argument and the explicit degree statement.

Let $\mathcal P_{\le m}(E)$ be polynomials of degree at most $m$.
Frobenius reciprocity gives

\[
 \operatorname{Hom}_{S_k}\bigl(V_\lambda,
             \mathbb R^k\otimes\mathcal P_{\le m}(E)\bigr)
 \simeq
 \operatorname{Hom}_{S_{k-1}}\bigl(
       \operatorname{Res}V_\lambda,\mathcal P_{\le m}(E)\bigr).
\]

The subgroup fixes the first label. On $E$, the other $k-1$
coordinates are independent, so its polynomial action is precisely
the permutation action in $k-1$ variables. The branching theorem
removes one corner box from $\lambda$, with multiplicity one:
Vershik--Okounkov, retained version math/0503040v3, Theorem 5.8,
printed page 22. Removing a box from row $i$ changes $n$ by $i-1$.
The bottom row is removable, hence the smallest degree on the right
is exactly $n(\lambda)-(r-1)$.

The degree-$m$ Hermite coefficient map in (S11) takes values in
$E\otimes\mathcal P_{\le m}(E)$, a submodule of the left-hand
target. It consequently vanishes below $D_\lambda$. Extra ambient
coordinates are fixed by $S_k$; their polynomial factors cannot
decrease this minimum active degree. Nonnegative Hermite coefficients
then prove the comparison in (S80). The argument bounds the entire
multiplicity space. It does not assert that the interface coefficient
map attains every degree allowed in the polynomial target.

### A quantitative gap for all sufficiently deep sectors

Write $\phi_0=\phi_1(0)$, $\sigma_\rho=\sqrt{1-\rho^2}$ and
$\ell=kv_k$. For every $\lambda$ with $D=D_\lambda\ge1$,
every $\rho\le b<1$ and every speed in that isotypic space,

\[
 \boxed{\displaystyle
 \mathfrak q_\rho(f)\ge
 \left[\frac{\ell}{\sigma_\rho}
  -\frac{4\phi_0}{b\sqrt{1-b^2}}
                    (\rho/b)^{D-1}\right]\|f\|^2.}
 \tag{S81}
\]

This gives a positive gap whenever its displayed bracket is positive;
the constant is independent of the ambient dimension. It does not
require a sign assertion for individual adjacent-facet terms.

To prove it, the absolute kernel row in (S10), before the $1/b$
factor, is
$2R_0+2\sum_lR_l\le2(2R_0+\sum_lR_l)=2a_b/b$,
using the standard row identity (S49). Gaussian differentiation of
$T_b(1_{T_1}-1_{T_2})$ in the unit normal direction gives
$0<a_b\le2b\phi_0/\sqrt{1-b^2}$: the derivative is bounded by
$(b/\sqrt{1-b^2})\mathbb E|Z|$. The symmetric Schur bound therefore
gives
$B_b(f)\le4\phi_0\|f\|^2/(b\sqrt{1-b^2})$.
Combine this with (S80) and the local lower bound
$a_\rho/\rho\ge\ell/\sigma_\rho$.

The constants can be made independent of Gaussian facet integrals.
With $Z_1,\ldots,Z_k$ independent standard normals,

\[
 \ell=\frac{\sqrt2}{k-1}\mathbb E\max_iZ_i
                  \ge\frac{2\phi_0}{k-1}.
\]

Indeed integrate out the common trace coordinate of the supporting
hyperplane to obtain
$v_k=\sqrt2\int\phi(t)^2\Phi(t)^{k-2}dt$.
Integration by parts in the maximum density
$k\phi(t)\Phi(t)^{k-1}$ gives the equality. The maximum increases
with $k$ and $\mathbb E\max(Z_1,Z_2)=1/\sqrt\pi$, giving the bound.
Thus the relative error in (S81), with $b=(1+\rho)/2$, is at most

\[
 \frac{8(k-1)}{\sqrt{(1+\rho)(3+\rho)}}
       \left(\frac{2\rho}{1+\rho}\right)^{D-1}
 <5(k-1)\left(\frac{2\rho}{1+\rho}\right)^{D-1}.
\]

In particular the gap is at least $\ell/(2\sigma_\rho)$ whenever

\[
 D\ge 1+\left\lceil
 \frac{\log(10(k-1))}{\log((1+\rho)/(2\rho))}\right\rceil.
\]

At fixed correlation, any sectors not covered by this bound have
$D_\lambda=O_\rho(\log k)$. If $s=r-1$ and
$t=k-\lambda_1$, then $D_\lambda\ge s(s-1)/2$ and
$t\le D_\lambda+s$. Hence these remaining partitions have
$O_\rho(\log k)$ boxes outside their first row. This is a restriction
on what remains to estimate, not a solution of those remaining blocks
or of the standard mass constraint.

A second convenient choice is $b=\sqrt{D/(D+1)}$, allowed when
$\rho\le b$. The square of the relative error is at most

\[
 4(k-1)^2(1-\rho^2)
       \frac{(D+1)^{D+1}}{D^D}\rho^{2D-2}.
\]

For example, an infinite family of previously untreated sectors now
has an explicit gap: for $k=2r$, $r\ge4$, $\lambda=(2^r)$ and
$0<\rho\le1/2$,

\[
 \mathfrak q_\rho(f)\ge
             \frac{2\ell}{3\sigma_\rho}\|f\|^2.
\]

Here $D=(r-1)^2\ge9$. The last squared error is increasing up to
$\rho=1/2$ on this interval. Since $(1+1/D)^D<3$, it is bounded by
$9(2r-1)^2((r-1)^2+1)/4^{(r-1)^2-1}$.
At $r=4$ this is $4410/65536<1/9$; on increasing $r$ its numerator
grows by a factor less than $4$ whereas its denominator grows by
$4^{2r-1}$. It therefore stays below $1/9$. The auxiliary estimate
$(1+1/D)^D<3$ follows by the binomial expansion and
$\sum_{j\ge0}1/j!<3$. This proof covers all $r\ge4$, all speeds
in the block and all ambient dimensions, rather than a finite scan
of cells or polynomials.

The checker independently computes graded character multiplicities
through $k=9$ to control the floor formula, checks the rational
constants and tests the first two-column instance. These are finite
algebraic controls. The unrestricted branching, Schur bound, Gaussian
maximum identity and spectral conclusions are the written proof.

### An explicit first-moment gap from convex covariance

The same degree structure gives an independent proof of the
first-moment interface theorem used in (S8). In the active space,
for every $k\ge3$ and balanced speed,

\[
 \boxed{Q_0(f)=\ell\|f\|^2-\|Df\|_F^2
       \ge\frac{\ell}{k}
                  \operatorname{dist}(f,\mathcal R_k)^2.}
 \tag{S82}
\]

Equality $Q_0(f)=0$ holds precisely for geometric rotation speeds.
This proof uses convex Gaussian covariance, the decomposition (S7)
and the Gaussian rotational identity (S55), rather than the accepted
global first-moment optimality theorem. In particular it does not
assume (S8) or its rotation gauge (S9).

Let $\mathsf B_0=D^*D$ on the balanced space. The degree-one target
in (S7) contains each nontrivial type at most once. Thus on each of
the standard, two-row and rotation-hook multiplicity spaces,
$\mathsf B_0$ has rank at most one. It vanishes on every other type,
also directly by (S80). The trivial target in (S7) has no interface
source. We now compute or bound these three possible eigenvalues.

On the balanced standard space, (S77) gives the rank-one eigenvalue
$\ell s_k^2/2\le\ell/2$. For the two-row type, set
$H^0=\{z\in E:z_1=z_2=0\}$ and $v_0(y)=P_0y$, with $P_0$ the
orthogonal projection onto $H^0$. The conditional facet law has

\[
 \mathbb E v_0=0,\qquad
 \mathbb E[v_0v_0^T]=\gamma_kP_0,\qquad0<\gamma_k\le1.
\]

The scalar covariance follows from the irreducible remaining-label
action; for $k=4$ it is a one-dimensional assertion. Positivity
follows from the open facet interior. The upper bound is the convex
Gaussian Poincare argument of (S77), applied to any unit direction
in $H^0$, with the same bounded smoothing and moment limit.

The two-row kernel is the scalar kernel (S57) on these vector fields.
Its first derivative uses

\[
 A=2I+\sum_{l=3}^kD_l
     =kI-\sum_lw_lw_l^T,\qquad A|_{H^0}=(k-1)I.
\]

Since $\int v_0d\nu=0$ and
$\int yv_0(y)^Td\nu=v_k\gamma_kP_0$, differentiation at zero
gives
$\mathsf B_0v_0=(k-1)v_k\gamma_kv_0$ in the complete two-row
multiplicity model. This nonzero field therefore identifies its
only nonzero multiplicity eigenvalue. It is at most
$(k-1)v_k=(1-1/k)\ell$, so this whole sector has first-moment
gap at least $\ell/k$. For $k=3$ this sector is absent.

Finally the rotation identity gives
$\mathsf B_0r=\ell r$ for $r\in\mathcal R_k$.
One may take the limit of the exact Gaussian identity (S55): on
balanced speeds the interaction tends in Hilbert--Schmidt norm to
$D^*D$, while $a_\rho/\rho$ tends strongly to $\ell$.
Both convergence assertions follow from the Gaussian Taylor
domination used in (S77), applied to the finitely many facets;
no positivity of $Q_0$ is used in this limit. The injective geometric
rotation map supplies the one rotation-hook copy. Since the target
has multiplicity one, $\mathsf B_0$ on that entire isotypic space
is exactly $\ell$ times the orthogonal projection onto this copy.
Every other copy consequently has gap $\ell$.

Combining the three bounds with orthogonality of representation
types proves (S82) and the asserted kernel. The rotation-gauge
identity (S8) now follows by polarization, with no global
first-moment optimality assumption. This is a local interface
statement; it does not prove global first-moment optimality over
arbitrary partitions.

### An explicit initial full gap without global first-moment optimality

For every $k\ge3$, every $d\ge k-1$, every balanced Gaussian $L^2$
facet speed and every $0<\rho\le[6\sqrt3\,k(k-1)]^{-1}$,

\[
 \boxed{\mathfrak q_{\rho,d}(f)\ge\frac{\ell}{2k}
                 \operatorname{dist}(f,\mathcal R_{k,d})^2.}
 \tag{S83}
\]

The interaction remainder has a direct bound, so no unevaluated
continuity modulus is needed. In the balanced compression the
constant Hermite term vanishes. Nonnegativity of every coefficient
in (S11) gives, for $0<\rho\le b<1$,

\[
 0\le\mathsf B_\rho-\mathsf B_0
       \le\frac\rho b(\mathsf B_b-\mathsf B_0)
       \le\frac\rho b\mathsf B_b,
 \qquad \mathsf B_0=D^*D.
\]

The absolute-row estimate in (S81), with $b=\sqrt{2/3}$, yields

\[
 \|\mathsf B_\rho-\mathsf B_0\|
        \le\frac{4\phi_0\rho}{b^2\sqrt{1-b^2}}
        =6\sqrt3\phi_0\rho\qquad(\rho\le\sqrt{2/3}).
\]

For balanced $g\perp\mathcal R_k$, combine $L_\rho\ge\ell I$
with (S82) to obtain

\[
 \mathfrak q_\rho(g)\ge
        (\ell/k-6\sqrt3\phi_0\rho)\|g\|^2.
\]

Since $\ell\ge2\phi_0/(k-1)$, the specified endpoint makes the
subtracted error at most $\ell/(2k)$. Exact rotational mixed
nullity (S55) allows subtraction of the rotation projection for
every $f$. The lift (S42) preserves this constant because
$\ell(1-\rho)/\sqrt{1-\rho^2}\ge\ell/\sqrt3>\ell/(2k)$
for $\rho\le1/2$. This proves (S83).

Both the gap and the interval endpoint are explicit. The multiplier
does not converge uniformly on the unbounded facets; only its global
lower bound is used. This route removes the accepted global
first-moment theorem from the initial interface theorem, without
establishing any new global optimality claim or an arbitrary-correlation
full gap. In particular (S83) is an explicit, restricted-correlation
extension of (S37) to every cell count.

## An exact two-variable kernel for tied-height standard speeds

A useful test family for the remaining ordinary-mass inequality has
an exact reduction independent of the number of coordinates. Realize
the cylindrical simplex in $\mathbb R^k$ with cells
$T_i=\{x_i\ge x_j\ \forall j\}$. The common-coordinate direction
is unused by the partition. On $F_{ij}$ take
$f_{ij}(x)=(c_i-c_j)h(x_i)$, where $c\in E\setminus\{0\}$.
The tied height is the actual coordinate $x_i=x_j$ in this
$k$-dimensional realization, rather than its projection onto $E$.

Let $\phi,\Phi$ be the standard normal density and distribution
function, and let $p_\rho(u,v),F_\rho(u,v)$ be the density and lower
distribution function of a standard bivariate normal pair with
correlation $\rho$. Write $F_u=\partial_uF_\rho$ and
$F_v=\partial_vF_\rho$, and define

\[
 m_k(u)=\sqrt2\,\phi(u)^2\Phi(u)^{k-2},\qquad
 \mathcal J_{k,\rho}(u,v)=
 4p_\rho(u,v)^2F_\rho(u,v)^{k-2}
 +2(k-2)p_\rho(u,v)F_u(u,v)F_v(u,v)F_\rho(u,v)^{k-3}.
\]

For $k\ge3$, $0<\rho<1$ and $h\in L^2(m_kdu)$, the exact norm,
balance constraint and form are

\[
 \boxed{\begin{aligned}
 \|f\|^2&=k|c|^2\int m_kh^2,\qquad
 Cf=0\ \Longleftrightarrow\ \int m_kh=0,\\
 \frac{\mathfrak q_\rho(f)}{k|c|^2}
   &=\iint\mathcal J_{k,\rho}(u,v)
            \bigl[h(u)^2-\rho^{-1}h(u)h(v)\bigr]du\,dv.
 \end{aligned}}
 \tag{S84}
\]

The second term in $\mathcal J$ is indispensable: it is the full
adjacent-facet contribution, not a tail correction to the same-facet
kernel. The formula includes every dimension of this height family
but does not cover every $G_0$-invariant function in (S51).

For a direct derivation, let $(X_l,Y_l)$ be independent standard
normal pairs of correlation $\rho$, $1\le l\le k$. Surface
integration on $X_1=X_2$ introduces
$\delta((X_1-X_2)/\sqrt2)=\sqrt2\delta(X_1-X_2)$.
Integrating the other $X$ coordinates below their common value $u$
gives $m_k(u)du$. Thus $\int m_k=v_k$.
For the same-facet Mehler interaction impose
$X_1=X_2=u$, $Y_1=Y_2=v$. The two surface factors multiply to $2$.
The two tied pairs supply $p_\rho^2$, and the other $k-2$ pairs
supply $F_\rho^{k-2}$. Multiplying by the coefficient $2$ in
(S49) gives the first term of $\mathcal J$.

For one adjacent term impose $X_1=X_2=u$, $Y_1=Y_3=v$.
The pair with label $1$ supplies $p_\rho(u,v)$, while integration
of $Y_2\le v$ and $X_3\le u$ supplies respectively

\[
 F_u=\phi(u)\Phi((v-\rho u)/\sigma),\qquad
 F_v=\phi(v)\Phi((u-\rho v)/\sigma).
\]

The other $k-3$ pairs supply $F_\rho^{k-3}$. The surface factor
is again $2$, and (S49) has $k-2$ such adjacent terms. This gives
the second term of $\mathcal J$. Integrating the row identity
$K_{\rm std}1=a_\rho/\rho$ against $h(u)^2$ supplies the local
term in (S84). The norm and balance follow from (S49) or directly
from $\sum_{i<j}(c_i-c_j)^2=k|c|^2$ and $\sum_ic_i=0$.
The argument first applies to bounded $h$; boundedness of the full
interface form and the displayed norm extend it to the stated
$L^2$ space. No numerical integration enters the reduction.

This is a one-dimensional operator problem with a two-variable kernel.
Set $r(u)=\int\mathcal J_{k,\rho}(u,v)dv$.
The local bounds give
$\ell/\sigma\le r(u)/m_k(u)\le2\phi_0/\sigma$.
On $L^2(du)$ the normalized interaction and the ordinary balance
vector are

\[
 T_h(u,v)=\frac{\mathcal J_{k,\rho}(u,v)}
                      {\rho\sqrt{r(u)r(v)}},\qquad
 b_h(u)=\frac{m_k(u)}{\sqrt{r(u)}}.
\]

The positive translation eigenvector is $\sqrt r$, with eigenvalue
$1/\rho$. The operator is compact, as a bounded compression and
normalization of the full compact interaction. Consequently the
two tests (S70) apply to this family, with $T_h,b_h$; the balance
condition is still ordinary $m_k$-mean, not $r$-mean. Positivity of
this restricted family would not prove (S51). A negative member,
however, would refute the full cylindrical interface assertion,
and therefore its active-space version by (S42).

Two identities check the normalization and connect this reduction
to the simplex noise stability $S_\rho(T)$:

\[
 \mathcal J_{k,0}(u,v)=k\,m_k(u)m_k(v),\qquad
 \iint\mathcal J_{k,\rho}=
                    \frac{2}{k(k-1)}\frac{d}{d\rho}S_\rho(T).
\]

For the second identity,
$S_\rho(T)=k\int p_\rho F_\rho^{k-1}$.
Direct differentiation of the Gaussian density gives
$\partial_\rho p_\rho=\partial_u\partial_vp_\rho$ and hence
$\partial_\rho F_\rho=p_\rho$. Integrating by parts twice shows
that the derivative of $\int p_\rho F_\rho^{k-1}$ is
$(k-1)/2$ times $\int\mathcal J$. Gaussian decay justifies the
differentiations and boundary limits at every fixed $0<\rho<1$.

The linear test $h(u)=u-\mu_k$, with
$\mu_k=v_k^{-1}\int um_k(u)du$, has an additional exact feature.
Write $x=x_E+z\mathbf1/\sqrt k$. Then
$h(x_i)=(x_{E,i}-\mu_k)+z/\sqrt k$.
The second term is a degree-one unused-coordinate standard speed,
which is a geometric rotation. Different unused Hermite degrees
are orthogonal for the form. Therefore (S84) for this linear test
equals the active-space standard form of $x_{E,i}-\mu_k$;
its extra trace contribution has exactly zero form. Nonlinear
height functions generally contain other unused-coordinate degrees,
so this identification is not asserted for them.

The checker controls the Gaussian conditioning matrices, surface
coefficients and formal differentiation coefficients in a stated
finite range. The all-$k$ surface integration, compact compression
and trace decomposition above are written proofs. Uncertified
height-family quadrature is not part of the proof receipt and
supplies no sign conclusion for (S84) or (S51).

## Triple antisymmetry isolates the remaining transport obstruction

The adjacent transports admit a simpler description that does not
depend on the shape of the representation. Let $R_\lambda$ be an
orthogonal irreducible representation, $P=(I-R_\lambda((12)))/2$
and $V^-=\operatorname{ran}P$. For $l\ge3$ put

\[
 A_l=\frac16\sum_{g\in S_{\{1,2,l\}}}
                       \operatorname{sgn}(g)R_\lambda(g).
\]

This is the orthogonal projection onto the subspace transforming
by the sign character of these three labels. Its range lies in
$V^-$. The exact full-facet multiplicity kernel is

\[
 \boxed{\begin{aligned}
 T_l&=2PR_\lambda((2l))P|_{V^-}=I-3A_l|_{V^-},\\
 K_{\lambda,\rho}(x,y)
   &=K_{{\rm std},\rho}(x,y)I
                         -3\sum_{l=3}^kM_l(x,y)A_l|_{V^-}.
 \end{aligned}}
 \tag{S85}
\]

Here the fields are $G_0$-equivariant on the full reference facet,
and $M_l=\Gamma_\rho e^{\tau x\cdot(2l)y}$ as before. To check
the algebra without an imported character calculation, write
$s=(12)$, $t=(2l)$ and $u=sts=(1l)$ in $\mathbb R[S_3]$.
Then
$2PtP=(t-st-ts+u)/2=P-3A_l$.
The self coefficient is $2I$, and the two adjacent facets give
$R(t)+R(sts)=2PR(t)P$ on $V^-$, which proves (S85).

For the rotation hook $A_l=(2/3)w_lw_l^T$, recovering $C_l$.
For the sign representation $A_l=I$, recovering the self-minus-
adjacent scalar kernel. Thus the possible negative adjacent
transport is exactly triple-label antisymmetry. In particular all
$A_l$ vanish if $\lambda$ has at most two rows. Repeated one-box
branching, in the retained Vershik--Okounkov Theorem 5.8, proves
this: restriction to three labels cannot contain the three-row
shape $(1,1,1)$. Conversely a shape with at least three rows admits
a Young-graph path beginning with this vertical three-strip, so
its three-label sign projection is nonzero. The same assertion
holds for any chosen three labels by conjugation.

The positive-semidefinite projectors in (S85) do not supply an
operator ordering of the integral kernels. For example on the ray
$x=sw_l$, $y=tw_l$ inside $F_{12}$, the exponential part of $M_l$
is $e^{-3\tau st/4}$. At two distinct positive ray parameters its
two-by-two kernel matrix has negative determinant. The common
Gaussian factors are positive row/column congruences, so they do
not change this sign. Continuous small supported test functions
give the corresponding failure of positive semidefiniteness on
the facet. Thus simply discarding the last term of (S85) is not
justified. The vanishing of that term in two-row representations,
however, can be used exactly.

## A single untied-coordinate reflection controls every nontrivial scalar type

Fix $k\ge4$, $0<\rho<1$ and $s=(34)$. Write
$F^+=\{x\in F_{12}:x_3>x_4\}$. The scalar standard kernel
commutes with $G_0=S_{\{3,\ldots,k\}}$. On scalar functions odd
under $s$, its folded kernel is

\[
 K^{\rm odd}(x,y)=K_{\rm std}(x,y)-K_{\rm std}(x,sy)
                                      \quad(x,y\in F^+).
\]

It is strictly positive on the relative interior of this product.
No ordering of the other coordinates is needed. The self term is
strictly positive, since
$x\cdot y-x\cdot sy=(x_3-x_4)(y_3-y_4)>0$.
Every adjacent term except those indexed by $3,4$ has the same
positive exponent difference. For the two exceptional terms, put
$E_l(y)=e^{\tau x\cdot(2l)y}$ and

\[
 \alpha=\tau(x_2-x_3)\ge0,\quad
 \beta=\tau(x_3-x_4)>0,\quad
 u=y_2-y_3\ge0,\quad t=y_3-y_4>0.
\]

The exact factorization already used in (S58) is

\[
 \frac{E_3(y)+E_4(y)-E_3(sy)-E_4(sy)}{E_3(sy)}
  =(e^{\alpha t}-1)(e^{\beta t}-e^{-\beta u})
                                +(e^{\beta t}-1)>0.
\]

This proves the stated positivity on the half-facet, rather than
only on the fully ordered chamber.

Set $w(x)=x_3-x_4>0$ on $F^+$. Projecting the Gaussian rotation
identity $K_1y=a_\rho x$ onto the coordinate difference $3-4$
gives exactly the component of (S59)

\[
 a_\rho w-K_{\rm std}w=2(R_3-R_4),\qquad
 R_l(x)=\int_{F_{12}}M_l(x,y)(y_2-y_l)d\nu_{12}(y).
\]

Here $R_3>R_4$ on all of $F^+$: substitute $y\mapsto sy$ in
$R_4$. Both integrals then have the factor $y_2-y_3$, and the
exponent difference is
$\tau(x_3-x_4)(y_2-y_4)>0$ in the facet interior. Since $w$ is
odd, $K_{\rm std}w$ is precisely the folded integral of
$K^{\rm odd}w$. Thus $\psi=\sqrt{a_\rho}w>0$ is a strict
supersolution for the compact self-adjoint positive-kernel operator
$T^{\rm odd}=a_\rho^{-1/2}K^{\rm odd}a_\rho^{-1/2}$:
$T^{\rm odd}\psi<\psi$ almost everywhere. The eigenvector
argument of (S73), pairing $|Tz|\le T|z|$ with $\psi$, gives
$|\lambda|<1$ for every nonzero eigenvalue. Compactness yields
$\theta_{k,\rho}=\|T^{\rm odd}\|<1$.

This also controls the entire nontrivial $G_0$ part of the scalar
space, not just one chosen odd subspace. Let $E_{G_0}$ denote group
averaging on $L^2(F_{12},\nu_{12})$. The coefficient $a_\rho$ is
$G_0$-invariant, so normalization commutes with the group action.
Every nontrivial irreducible $G_0$ representation has a nonzero
$s$-minus space: if one transposition acted as the identity on the
whole representation, conjugacy and generation would make the
representation trivial. Every nonzero eigenvalue of the compact
normalized scalar operator on $\ker E_{G_0}$ consequently also
occurs on the $s$-odd space. Its norm is therefore at most
$\theta_{k,\rho}$. It follows that

\[
 \boxed{\frac1\rho\left[\int a_\rho h^2d\nu_{12}
                  -\iint K_{\rm std}hh\,d\nu_{12}^2\right]
 \ge\frac{\ell}{\sigma}(1-\theta_{k,\rho})\|h\|^2
      \quad\text{if }E_{G_0}h=0.}
 \tag{S86}
\]

The remaining scalar standard problem (S51) instead has
$E_{G_0}h=h$ and only its ordinary Gaussian mean zero. It is not
covered by (S86). The invariant translation function $h=1$ has
normalized eigenvalue $1/\rho>1$, so removal of the nontrivial-type
condition would be false.
In contrast $E_{G_0}h=0$ implies both $\int h\,d\nu_{12}=0$ and
$\int a_\rho h\,d\nu_{12}=0$, because the measure and $a_\rho$
are group invariant. The mass-weight mismatch in (S70) is therefore
absent on the scalar types controlled by this reflection.

## Every nonstandard two-row block has an all-correlation gap

For every $k\ge4$, $2\le p\le\lfloor k/2\rfloor$,
$d\ge k-1$ and $0<\rho<1$, the complete
$V_{(k-p,p)}$ isotypic interface space has no kernel and satisfies

\[
 \boxed{\mathfrak q_{\rho,d}(f)
       \ge\kappa^{\rm 2row}_{k,\rho}\|f\|^2,\qquad
 \kappa^{\rm 2row}_{k,\rho}=
                \frac{\ell}{\sigma}(1-\theta_{k,\rho})>0.}
 \tag{S87}
\]

The constant is independent of both the ambient dimension and the
choice of $p$; its spectral factor is qualitative. This includes
all multiplicities and all $L^2$ speeds. It extends (S60), whose
case is $p=2$, and includes $(3,3)$ at six cells and every higher
nonstandard two-row shape.

Indeed the vertical two-strip rule (S4) identifies $V^-_{(k-p,p)}$
under $G_0$ with the single irreducible
$V_{(k-p-1,p-1)}$. Removing one box from each of the two rows is
the only allowed removal. Formula (S85) makes the raw kernel
$K_{\rm std}I$ on this entire multiplicity space. For $p\ge2$
the remaining-label representation is nontrivial. Every component
of an equivariant field therefore has zero group average, so
(S86), applied componentwise, proves the active-space assertion.
The constant normalization of the multiplicity norm cancels.
For $p=1$ the remaining-label representation is trivial, exactly
the standard sector excluded from this argument.

For extra unused coordinates, orthogonal Hermite degree $m$
multiplies the normalized interaction by $\rho^m$. Its norm is at
most $\rho^m\theta_{k,\rho}\le\theta_{k,\rho}$, with unchanged
local coefficient. Parseval gives (S87) in every ambient dimension.
No mass constraint is present in these nonstandard types.

The checker verifies the three-label group-algebra corner and its
Young-model projections in a finite range, the exponential
factorization and selected half-facet comparisons without ordering
the other coordinates. It rejects the wrong half-facet and the
positive-projector/operator-ordering shortcut. The unrestricted
branching, strict supersolution, nontrivial-type spectral transfer
and whole-block theorem are written internal proofs, without
independent replay, external review or formal verification.

## Antisymmetry gives an exact hierarchy by representation height

For a label set $S$, write
$E_S=(|S|!)^{-1}\sum_{g\in S_S}\operatorname{sgn}(g)R_\lambda(g)$.
These are orthogonal sign projections, with $E_{\{1,2\}}=P$ and
$E_{\{1,2,l\}}=A_l$. If $|S|=m$, $r\notin S$ and
$t=(jr)$ for any $j\in S$, the exact group-algebra identities are

\[
 E_SR(t)E_S=\frac1mE_S-\frac{m+1}{m}E_{S\cup\{r\}},\qquad
 E_SE_{tS}E_S=\frac1{m^2}E_S+
                 \left(1-\frac1{m^2}\right)E_{S\cup\{r\}}.
\]

Indeed right cosets of $S_S$ in $S_{S\cup\{r\}}$ give
$(m+1)E_{S\cup\{r\}}=E_S-\sum_{i\in S}E_SR((ir))$.
Multiplying on the right by $E_S$, all $m$ compressed transpositions
are equal by conjugation inside $S_S$. This proves the first identity.
Squaring it gives the second because the larger sign projection is
nested in $E_S$. For three-label projections in a representation with
at most three rows this gives $A_lA_rA_l=A_l/9$ for distinct $l,r$.
The possible four-label sign subspace is exactly the correction to
that relation.

There is also an exact Gaussian-kernel hierarchy. Take
$S=\{1,\ldots,m\}$ with $2\le m\le k$, and put $u_S=E_Sh$ for
a full-facet $G_0$-equivariant $V^-$-valued field $h$. Define the
scalar kernel

\[
 L_m(x,y)=2M_0(x,y)-2\sum_{r=3}^mM_r(x,y)
  +\sum_{l=m+1}^k\left[M_l(x,y)
                  -\sum_{r=3}^mM_l(x,(rl)y)\right].
\]

Then

\[
 \boxed{E_SK_\lambda h
   =L_m(E_Sh)-(m+1)\sum_{l=m+1}^k
       \int_{F_{12}}M_l(\,\cdot\,,y)
                        E_{S\cup\{l\}}h(y)d\nu_{12}(y).}
 \tag{S88}
\]

To verify every term, for $l\in S\setminus\{1,2\}$ the transport
has $E_ST_l=-2E_S$. For $l\notin S$, multiply the preceding coset
identity by $P$. The two terms indexed by $1,2$ combine to $E_ST_l$,
and all the other transpositions belong to $G_0$. Thus

\[
 E_ST_lh(y)=u_S(y)-\sum_{r=3}^mu_S((rl)y)
                         -(m+1)E_{S\cup\{l\}}h(y).
\]

Changing $y\mapsto(rl)y$ in those integrals proves (S88), with
unchanged surface measure and facet. It first holds for bounded
fields and extends by the bounded Gaussian operators. At $m=2$
it is (S85); at $m=k$ it is the full sign kernel. If $\lambda$ has
exactly $m$ rows, every higher sign projection in (S88) vanishes
by branching, and the projected field satisfies a scalar intertwining.
This statement by itself does not assert that the highest projection
detects every field in an arbitrary-height multiplicity space.

## All three-row blocks embed into one scalar spectral problem

The triple projections do detect every higher-row field, with a
quantitative bound. Put $q=k-2$ and
$c_\lambda=\sum_{(i,j)\in\lambda}(j-i)$, the sum of box contents.
On a remaining-label branch $V_\mu\subset V_\lambda^-$ from (S4),

\[
 \sum_{l=3}^kA_l
       =\frac{k-3-c_\lambda+c_\mu}{3}I.
\]

For completeness, let $C_k$ be the sum of all transpositions and
$C_{G_0}$ the sum of those on the remaining labels. On $V^-$,
$C_k=-I+\sum_lT_l+C_{G_0}$; substitute (S85).
The classical content eigenvalues of the transposition sums give
the displayed formula. They follow by summing the Young--Jucys--Murphy
elements: Vershik--Okounkov, retained version math/0503040v3,
definition and (2.1) on printed page 10, and the content-spectrum
statement in Theorem 5.8 on page 22. These exact statements were
inspected; no Gaussian spectral comparison is imported from them.

If $\lambda/\mu$ removes boxes from rows $i<j$, this eigenvalue is

\[
 \frac{\sum_{r\notin\{i,j\}}\lambda_r+i+j-3}{3}.
\]

It is zero in two-row shapes and is at least $1/3$ in every shape
with three or more rows. Hence these sign subspaces span $V^-$.
For equivariant fields in any such higher-row block this also gives
$\|A_3h\|^2\ge\|h\|^2/(3q)$ by group averaging.
For $\lambda=(a,b,c)$ with three positive rows, the possible
eigenvalues are $c/3,(b+1)/3,(a+2)/3$, omitting any inadmissible
strip removal. In particular, for $G_0$-equivariant fields,

\[
 \frac{c}{3q}\|h\|^2\le\|A_3h\|^2
        =\frac1q\sum_{l=3}^k\|A_lh\|^2
                     \le\frac{a+2}{3q}\|h\|^2.
\]

The middle equality uses $A_l=R((3l))A_3R((3l))$ and a
measure-preserving change of variables. Thus $h\mapsto A_3h$ is
a bounded injection with closed image on the full multiplicity space.
The same bounds hold after normalization by $\sqrt{a_\rho}$,
since this scalar multiplier is $G_0$-invariant.

Set $\mathcal L_3=a_\rho^{-1/2}L_3a_\rho^{-1/2}$ on scalar
$L^2(F_{12},\nu_{12})$. The explicit kernel is

\[
 L_3=2M_0-2M_3+
                \sum_{l=4}^k\bigl[M_l-M_l(\,\cdot\,,(3l)\,\cdot\,)\bigr].
\]

It is compact but is not asserted to be self-adjoint in this scalar
inner product. For every three-row $\lambda$, (S88) and the closed
injection give

\[
 \boxed{A_3(a_\rho^{-1/2}K_\lambda a_\rho^{-1/2})
                =\mathcal L_3 A_3.}
 \tag{S89}
\]

In particular every nonzero eigenvalue of every normalized three-row
block is an eigenvalue of this one scalar operator: a nonzero
component of $A_3h$ supplies a scalar eigenfunction. The positive
function $\sqrt{a_\rho}(x_2-x_3)$ has eigenvalue $1$, by applying
the intertwining to the exact rotation field.

The three-label sign fiber transforms under
$H_3=S_{\{4,\ldots,k\}}$ by the shape $(a-1,b-1,c-1)$, from the
vertical three-strip rule. Consequently the scalar components lie
in $H_3$ types with at most three rows. The scalar kernel commutes
with $H_3$. Let $\mathscr U_3$ be the closed sum of those isotypic
subspaces of scalar $L^2(F_{12},\nu_{12})$. It would suffice to
bound the spectral radius of $\mathcal L_3|_{\mathscr U_3}$ by $1$
and prove that its fixed space consists only of the positive
function above. The sign fiber has a fixed vector only for the
rotation hook. Thus these scalar assertions would exclude eigenvalue
$1$ in every other three-row block; compactness would give a strict
gap. Neither scalar assertion is proved here. A bound on the whole
scalar space is stronger than this sufficient condition, and no
corresponding sufficiency claim for every higher-row block follows
from (S89).

The positive-kernel argument of (S86) cannot be repeated for this
scalarization. Already at $k=5$, take

\[
 x=(7,7,-6,0,-8),\qquad y=(27,27,7,-23,-38)/5.
\]

Both lie in the facet interior and have $x_4>x_5$, $y_4>y_5$.
At $\tau=5\log2$, the difference
$[L_3(x,y)-L_3(x,(45)y)]/\Gamma_\rho(x,y)$ has leading term
$-2^{680}$. Its positive terms have coefficient sum $8$ and
exponents at most $650$, so it is less than
$-2^{680}+8\,2^{650}<0$.
Rescaling both points preserves this sign at every fixed positive
correlation. This is a pointwise obstruction to that particular
reflection proof, not an interface instability: the five-cell mixed
block is already positive by (S74). The spectral problem may require
a different cone or a direct operator estimate.

The exact checker verifies nested coset identities in specified finite
groups, the transposition-content frame identities and spectral
polynomials in finite Young models, and the displayed negative
reflection witness. The general coset proof, branching/content
application, bounded Gaussian intertwining and closed-image spectral
inclusion are the unrestricted written arguments. No new all-correlation
gap is inferred from this reduction.

### Highest antisymmetry detects every representation height

The three-row detection argument has an all-height extension. The
ingredients are the classical vertical/horizontal Pieri rule for every
strip size (Macdonald, Chapter I, (5.16)--(5.17), pp. 73--74, with the
characteristic map (7.3), (7.5)--(7.6), pp. 113--114) used in (S4),
and the Young branching and positive orthogonal-generator formulas
in the retained Vershik--Okounkov Theorem 5.8 and Section 6. The following
frame argument concerns these finite representations; none of those
imported statements supplies a Gaussian spectral estimate.

**Theorem S90 (highest-sign frame and scalar reduction).**

Let
$\lambda\vdash k$ have at least two rows, let $2\le m\le k$, and put

\[
 N_m=\binom{k-2}{m-2},\qquad
 F_{\lambda,m}=\left.\sum_{\substack{S\supset\{1,2\}\\|S|=m}}
                           E_S\right|_{V_\lambda^-}.
\]

If $\lambda$ has fewer than $m$ rows, $F_{\lambda,m}=0$.
If it has at least $m$ rows, this finite-dimensional self-adjoint
operator is strictly positive. Write
$\alpha_{\lambda,m}=\min\operatorname{spec}F_{\lambda,m}>0$.
Every $G_0$-equivariant $L^2$ field then satisfies

\[
 \boxed{\quad
 \|E_{[m]}h\|^2=\frac1{N_m}\int
          \langle h,F_{\lambda,m}h\rangle\,d\nu_{12}
       \ge\frac{\alpha_{\lambda,m}}{N_m}\|h\|^2.
 \quad} \tag{S90}
\]

Here $[m]=\{1,\ldots,m\}$. The same assertion holds with any
positive $G_0$-invariant scalar weight, in particular $a_\rho$.
The constant is qualitative and finite-dimensional. A single $E_{[m]}$
need not be injective on the pointwise vector space; (S90) uses
equivariance on the entire facet.

To prove strict positivity, $F_{\lambda,m}$ commutes with $G_0$.
By vertical two-strip Pieri,
$V_\lambda^-|_{G_0}$ is a multiplicity-free sum of $V_\mu$ with
$\lambda/\mu$ a vertical two-strip. It suffices to show that each
such summand has nonzero overlap with $\operatorname{ran}E_{[m]}$.
Extend its two removed boxes to a vertical $m$-strip
$\lambda/\gamma$. This extension always exists when
$m\le\operatorname{height}(\lambda)$: within each block of equal
row lengths, a removable vertical strip selects a bottom suffix;
increase those suffixes one row at a time until $m$ rows are selected.
Thus $\gamma\subseteq\mu$ and $\mu/\gamma$ is a vertical
$(m-2)$-strip.

Twist by the sign character, replacing $\lambda,\mu,\gamma$ by
their conjugates and $E_{[m]}$ by the $S_{[m]}$ symmetrizer.
Use the orthogonal Young basis with the reversed label order
$k,k-1,\ldots,1$. There is a basis tableau whose first $k-m$
boxes have shape $\gamma'$ and first $k-2$ boxes have shape
$\mu'$: complete these three nested diagrams successively.
It belongs to the $G_0$ isotypic space of shape $\mu'$.
Fix the initial filling of $\gamma'$. On the remaining skew-tableau
space, $\lambda'/\gamma'$ is a horizontal $m$-strip. Horizontal
Pieri says that its $S_{[m]}$ invariant space has dimension one.
The graph of its standard skew tableaux under admissible adjacent
swaps is connected: these tableaux are the linear extensions of the
box order, and moving the first differing element to its prescribed
place connects any two linear extensions by swaps of incomparable
neighbors.

In this orthogonal basis each $I+R(s_i)$ is entrywise nonnegative,
with strictly positive entries on those edges. The average of the
projections $(I+R(s_i))/2$, $1\le i<m$, is a nonnegative
irreducible matrix with spectrum in $[0,1]$. Its eigenvalue $1$ is
simple by horizontal Pieri, so its invariant vector is strictly
positive. Its powers tend to the symmetrizer, whose diagonal at our
chosen tableau is therefore positive. Consequently the trace of
the $G_0$ isotypic projection times $E_{[m]}$ is positive.
The latter projector has range in $V_\lambda^-$ before twisting.
Every conjugate $E_S$ has the same trace on that summand, proving
that $F_{\lambda,m}$ has a positive scalar there. This proves
strict positivity on all of $V_\lambda^-$.

Finally $G_0$ acts transitively on the $N_m$ sets $S$, and its
action preserves the facet measure and the scalar weight.
Equivariance makes all $\|E_Sh\|$ equal. Summing their squares
gives (S90), with no pointwise injectivity assumption.

Suppose now that $\lambda$ has **exactly $m$ rows**. The next
antisymmetrizers in (S88) vanish. Hence, with
$\mathcal L_m=a_\rho^{-1/2}L_ma_\rho^{-1/2}$,

\[
 E_{[m]}\bigl(a_\rho^{-1/2}K_\lambda a_\rho^{-1/2}\bigr)
       =(\mathcal L_m\otimes I)E_{[m]}.
\]

This is a closed-image similarity, not just a formal projection.
Its image can be specified without taking a closure. Put
$B_m=S_{\{3,\ldots,m\}}\times S_{\{m+1,\ldots,k\}}$ and
choose representatives $g$ of $G_0/B_m$ carrying
$\{3,\ldots,m\}$ to each $(m-2)$-subset. For a field $u$ with
values in $E_{[m]}V_\lambda$ and
$u(bx)=R_\lambda(b)u(x)$, define

\[
 (\mathcal R_{\lambda,m}u)(x)
   =F_{\lambda,m}^{-1}\sum_{g\in G_0/B_m}
                          R_\lambda(g)u(g^{-1}x),\qquad
 \Pi_{\lambda,m}=E_{[m]}\mathcal R_{\lambda,m}.
\]

The sum is independent of the representatives and is $G_0$-equivariant.
For $u=E_{[m]}h$ it equals $F_{\lambda,m}h$ before applying the
inverse, so $\mathcal R_{\lambda,m}E_{[m]}h=h$.
In fact $\Pi_{\lambda,m}$ is the orthogonal projection onto the
physical image: for $A:h\mapsto E_{[m]}h$, the adjoint is
$A^*u=N_m^{-1}\sum_gR_\lambda(g)u(g^{-1}x)$ and
$A^*A=F_{\lambda,m}/N_m$. Thus
$\Pi_{\lambda,m}=A(A^*A)^{-1}A^*$.
It commutes with multiplication by $a_\rho$.
The normalized block is similar, by the bounded bijection $A$, to
$(\mathcal L_m\otimes I)|_{\operatorname{ran}\Pi_{\lambda,m}}$.

The fiber $E_{[m]}V_\lambda$ is the sign representation of
$S_{\{3,\ldots,m\}}$ times the irreducible representation of
$H_m=S_{\{m+1,\ldots,k\}}$ with shape
$\nu=(\lambda_1-1,\ldots,\lambda_m-1)$, omitting zeros.
This is again vertical Pieri: a vertical strip of size $m$ removes
one box from every row. The scalar kernel $L_m$ commutes with both
factors of $B_m$. Consequently every nonzero eigenvalue of every
height-$m$ block occurs in this single scalar operator, on functions
alternating in labels $3,\ldots,m$ and in an $H_m$ type of at most
$m$ rows. Exact equivalence holds on the physical images above;
an estimate on the entire larger scalar space is only a sufficient,
possibly stronger condition.

At $m=2$ this recovers the two-row scalarization, and at $m=k$ it
recovers the sign block and its folded determinant. At $m=3$,
(S89) supplies sharper explicit frame constants. No scalar spectral
bound is proved by (S90). In particular the ordinary balance constraint
in the standard block remains (S51), and the other unsolved
higher-row bounds remain necessary. The exact checker tests conjugate
frames in specified finite Young models and compatible vertical-strip
and reverse-tableau chains through ten cells; the unrestricted frame,
reconstruction and Gaussian similarity arguments are written proofs.

### Positive harmonic fields for every higher exterior hook

There is a common polynomial field behind the linear supersolution
in (S73) and the sign Vandermonde in (S63). Its positivity is valid in
every exterior degree. Its Gaussian facet integrals reduce to one
dimension, including the adjacent transports. The remaining comparison
between those integrals is stated separately below.

**Theorem S91 (positive hook fields and exact moment reduction).**

Let $k\ge4$, $q=k-2$, $2\le p\le q$, and $0<\rho<1$, with
$\sigma=\sqrt{1-\rho^2}$. Put
$e=p(p-1)/2$, $s_p=(-1)^{(p-1)(p-2)/2}$ and let $P_H$ be
orthogonal projection onto $H$. Define

\[
 W_p(x)=s_p\,p_0\wedge P_Hx\wedge P_H(x^{\odot2})
                      \wedge\cdots\wedge P_H(x^{\odot(p-1)}).
 \tag{S91a}
\]

Powers here are coordinate powers. This is a $G_0$-equivariant,
homogeneous harmonic polynomial field of degree $e$ on $H$. Every
coordinate in the increasing Helmert exterior basis is strictly
positive on $W$. It is invariant under adding a common constant to
the $k$ coordinates before projecting the power columns.
The corresponding full symmetry block is
$V_{(k-p-1,1^{p+1})}$, with the complete field model from (S72).

Here is a direct proof of the strict sign. Put
$X_j^{(r)}=\sum_{i=1}^j x_i^r-jx_{j+1}^r$ and
$\alpha_j=1/\sqrt{j(j+1)}$. The coordinate for
$2\le j_1<\cdots<j_p\le k-1$ is
$s_p\prod_t\alpha_{j_t}$ times the determinant with rows
$(1,X_{j_t}^{(1)},\ldots,X_{j_t}^{(p-1)})$.
Subtract consecutive rows, starting at the bottom, and use

\[
 X_{j_{t+1}}^{(r)}-X_{j_t}^{(r)}
   =\sum_{l=j_t}^{j_{t+1}-1}(l+1)
              \int_{x_{l+2}}^{x_{l+1}}r z^{r-1}\,dz.
\]

Multilinearity expresses the remaining determinant as $(p-1)!$
times a positive weighted sum of integrals of
$\det[z_t^{r-1}]_{t,r=1}^{p-1}$ over successively lower, disjoint
intervals. Almost everywhere in each such product,
$z_1>\cdots>z_{p-1}$. That Vandermonde has sign $s_p$ and is
nonzero. This proves the asserted strict positivity for every $p$.

For the harmonic and integral statements, set
$b_i=e_i-k^{-1}\mathbf1\in H$, $3\le i\le k$. Since
$\sum_{i=3}^k b_i=-2p_0$, expansion of the wedge gives the identity

\[
 W_p(x)=-\frac{s_p}{2}
       \sum_{\substack{I\subset\{3,\ldots,k\}\\|I|=p}}
          \det[x_i^r]_{\substack{i\in I\ \mathrm{increasing}\\0\le r<p}}
                      \bigwedge_{i\in I\ \mathrm{increasing}}b_i.
 \tag{S91b}
\]

Each determinant is a Vandermonde, independent of a common coordinate
shift and harmonic in its variables. For example, its Laplacian is
alternating of degree $e-2$, whereas every alternating polynomial is
divisible by the degree-$e$ Vandermonde. These observations also prove
harmonicity on $H$: use independent coordinates with the tied top
coordinate of variance $1/2$ and the other coordinates of variance
$1$, then remove the independent common trace. Thus for an isotropic
standard Gaussian $Z_H$ on $H$,
$\mathbb E W_p(\rho x+\sigma Z_H)=\rho^e W_p(x)$.

The same expansion computes the truncated integral exactly. For a
mean vector $\mu\in H$ define

\[
 T_r(\mu_i,u)=\int_{-\infty}^u z^r\,
               \frac1\sigma\phi((z-\mu_i)/\sigma)\,dz,\qquad
 g_{\mu_2}(u)=\frac{\sqrt2}{\sigma}
                      \phi(\sqrt2(u-\mu_2)/\sigma).
\]

For a fixed linear map $C:H\to H$, write

\[
 \begin{split}
 P(\mu)&=\int_{\mathbb R}g_{\mu_2}(u)
                         \prod_{i=3}^kT_0(\mu_i,u)\,du,\\
 J_p(C,\mu)&=-\frac{s_p}{2}\sum_{|I|=p}
          \left(\bigwedge_{i\in I}Cb_i\right)
       \int_{\mathbb R}g_{\mu_2}(u)
          \det[T_r(\mu_i,u)]_{\substack{i\in I\\0\le r<p}}
                     \prod_{j\notin I,\ j\ge3}T_0(\mu_j,u)\,du.
 \end{split} \tag{S91c}
\]

All set indices in the wedge and determinant are increasing.
If $Y=\mu+\sigma Z_H$, these are respectively
$\mathbb P(Y\in F_{12})$ and
$\mathbb E[\mathbf1_{F_{12}}(Y)\Lambda^pC W_p(Y)]$.
Indeed, before removing the common trace, the tied coordinate
$U\sim N(\mu_2,\sigma^2/2)$ is independent of the coordinates
$Y_i\sim N(\mu_i,\sigma^2)$, $i\ge3$. Both the facet event and
(S91b) are invariant under that trace removal. Conditional on $U=u$,
the coordinates are independent. Each determinant term has distinct
coordinate indices, so its integral is the displayed determinant of
unnormalized moments, times the spectator probabilities. Polynomial
Gaussian integrability justifies the finite expansion and Fubini.
No multivariate quadrature is required for this identity.

The moments use only $\phi,\Phi$ and a finite recurrence. With
$z=(u-\mu_i)/\sigma$,

\[
 T_0=\Phi(z),\quad T_1=\mu_i\Phi(z)-\sigma\phi(z),\quad
 T_r=\mu_iT_{r-1}+(r-1)\sigma^2T_{r-2}
                              -\sigma u^{r-1}\phi(z)\quad(r\ge2).
\]

Let $\mu_0=\rho x$, $\mu_l=\rho D_lx$ and
$c_l=\exp[-\rho^2(w_l\cdot x)^2/(4\sigma^2)]$.
Completing the Gaussian square in each kernel gives the exact formulas

\[
 \begin{split}
 \mathsf K_pW_p(x)
   &=\frac{\phi_0}{\sigma}
          \left[2J_p(I,\mu_0)+\sum_{l=3}^kc_lJ_p(C_l,\mu_l)\right],\\
 a_\rho(x)/\rho
   &=\frac{\phi_0}{\sigma}
             \left[2P(\mu_0)+\sum_{l=3}^kc_lP(\mu_l)\right].
 \end{split} \tag{S91d}
\]

In particular the following is a concrete **sufficient comparison**
for every remaining exterior hook:

\[
 2J_p(I,\mu_0)+\sum_lc_lJ_p(C_l,\mu_l)
    \ \le\ \rho^e\left[2P(\mu_0)+\sum_lc_lP(\mu_l)\right]W_p(x)
 \tag{S91e}
\]

coordinatewise on $W$. The case $p=3$ is proved in (S92) below;
the intermediate cases $4\le p<q$ remain open. Its correlation
parameter can be removed exactly.
Write $P^{[1]},J_p^{[1]}$ for (S91c) at variance parameter $\sigma=1$,
and set $z=\rho x/\sigma$. Homogeneity and the substitution
$u=\sigma v$ give

\[
 P(\mu)=P^{[1]}(\mu/\sigma),\qquad
 J_p(C,\mu)=\sigma^eJ_p^{[1]}(C,\mu/\sigma).
\]

Consequently (S91e) is equivalent to the single unit-variance inequality

\[
 2J_p^{[1]}(I,z)+\sum_l e^{-(w_l\cdot z)^2/4}
                                      J_p^{[1]}(C_l,D_lz)
 \ \le\ \left[2P^{[1]}(z)+\sum_l e^{-(w_l\cdot z)^2/4}
                                      P^{[1]}(D_lz)\right]W_p(z).
 \tag{S91f}
\]

As $x$ ranges over $W$, so does $z$. Thus a proof of the **pointwise
field comparison for all $x$** at any one positive correlation would
prove it at every positive correlation. A spectral bound at one
correlation alone does not supply that hypothesis.
If (S91e) holds at a fixed $(k,p,\rho)$, the positive
folded kernel from (S72), applied to $\psi=\sqrt{a_\rho}W_p$,
gives $\|a_\rho^{-1/2}\mathsf K_pa_\rho^{-1/2}\|
\le\rho^{e-1}$. To see this, integrate
$|Tz|\le T|z|$ against $\psi$ for each eigenvector of the compact
self-adjoint $T$, and use (S91e). For $p\ge3$ this would yield the
explicit full-multiplicity gap
$\ell\sigma^{-1}(1-\rho^{e-1})$, independently of the ambient
dimension by the unused-coordinate Hermite factorization.
At $p=2$, (S73) already proves the strict residual behind the
non-strict factor $1$. At $p=q\ge3$, the comparison and explicit factor
are already proved by (S63)--(S66), and (S92) proves $p=3$.
The remaining intermediate-degree comparison is the obligation;
positivity and free harmonicity alone do not
control the truncated adjacent terms in (S91d).

The exact controls below check selected positive wedge coordinates,
(S91b), free Gaussian Vandermonde moments, the moment recurrence and
adjacent $C_lb_i$ identities. The all-degree positivity, harmonicity,
one-dimensional integration formula and conditional Schur conclusion
are written arguments. No finite numerical sampling certifies (S91e).

### The fourth exterior hook: a convex clipping argument

**Theorem S92 (fourth exterior hook).**

For every $k\ge5$, $d\ge k-1$ and $0<\rho<1$, the complete
$V_{(k-4,1^4)}$ isotypic block of the cylindrical interface form has
no kernel and satisfies

\[
 \boxed{\qquad \mathfrak q_\rho(f)\ge
             \ell\sqrt{1-\rho^2}\,\|f\|^2.\qquad} \tag{S92}
\]

This includes all $L^2$ multiplicities. At five cells it agrees with
the sign estimate; the new cases are $k\ge6$. It proves (S91e)--(S91f)
for $p=3$. The other intermediate exterior degrees and the standard
constraint are not resolved by this theorem.

#### Combining the conditional adjacent terms

Use unit variance as in (S91f), and subtract the common top mean.
Write the remaining means as $0>m_3>\cdots>m_k$, and put
$z=P_H(0,0,m_3,\ldots,m_k)$. Fix a common truncation level $u$.
Define

\[
 f(t)=\lambda(u-t),\quad \lambda(v)=\phi(v)/\Phi(v),\quad
 f_0=f(0),\quad f_i=f(m_i),\quad r_i=f_i/f_0,
\]

and let $M_{ir}$ be the $r$th moment of $N(m_i,1)$ conditioned to be
at most $u$, with $M_{0r}$ the corresponding moment at mean zero.
Set

\[
 A=2p_0+\sum_{i=3}^kr_i b_i,\qquad
 U_r=\sum_{i=3}^k(M_{ir}-M_{0r})b_i.
\]

For every $p\ge2$ the self and adjacent integrands in (S91f), after
division by the common positive factor
$g_0(u)\prod_i\Phi(u-m_i)$, sum exactly to

\[
                 s_p A\wedge U_1\wedge\cdots\wedge U_{p-1}.
 \tag{S92a}
\]

The scalar row integrand divided by that factor is $2+\sum_i r_i$.
Indeed the adjacent scalar density ratio is

\[
 e^{-m_l^2/4}\frac{g_{m_l/2}(u)\Phi(u)}
                         {g_0(u)\Phi(u-m_l)}=r_l.
\]

Under that adjacent law only coordinate $l$ changes from mean $m_l$
to mean zero. Since $C_lp_0=b_l$, $C_lb_l=-\sum_i b_i$ and
$C_lb_i=b_i$ for $i\ne l$, its conditional field expectation is
$s_p b_l\wedge U_1\wedge\cdots\wedge U_{p-1}$.
One can take expectations columnwise here: repeated occurrences of
one coordinate in the wedge have repeated vector factors and vanish;
the surviving terms have distinct independent coordinate indices.
Replacing $\sum_iM_{ir}b_i$ by $U_r$ in the self term adds only a
multiple of $p_0$. Summing proves (S92a). This identity does not assert
the higher-degree comparison by itself.

For $p=3$, let

\[
 F=P_H(f_0,f_0,f_3,\ldots,f_k),\quad
 Z=P_H(0,0,m_3^2,\ldots,m_k^2),\quad
 T=P_H(0,0,m_3f_3,\ldots,m_kf_k),\quad
 S=2f_0+\sum_i f_i.
\]

The first two truncated-moment formulas give
$A=F/f_0$, $U_1=z-F$ and $U_2=Z-T-uF$.
Consequently the desired conditional residual, multiplied by $f_0$,
is

\[
 R(f)=S W_3(z)+F\wedge z\wedge(Z-T)
       =L(f)+Q(f),\qquad
 L(f)=-(Sp_0-F)\wedge z\wedge Z,\quad
 Q(f)=-F\wedge z\wedge T. \tag{S92b}
\]

We prove that every Helmert coordinate of $R(f)$ is nonnegative for
any positive increasing convex function $f$ on $[m_k,0]$ whose
secant slopes lie in $[0,1]$.

#### A three-coordinate clipping inequality

For $0<a<b<c$ write $X=(a,b,c)^T$ and
$\det(v,X,w)$ for the determinant with those three columns.
Let $Y$ be a convex piecewise linear function on $[0,\infty)$ with
$Y(0)=0$. For a threshold $S\ge0$ write
$A_S=\min(X,S)$, $V_S=S\mathbf1-A_S$ and
$Y_S=Y(\min(X,S))$, all entrywise. Then

\[
 L_S=\det(V_S,X,Y(X))\ge0,\qquad
 D_{S,T}=\det(A_S,X,A_T)\le0\quad(0\le S\le T),
\]

and

\[
 C_0=(S+T)\det(\mathbf1,X,Y(X))
       -\det(A_S,X,Y(X)-Y_T)
       -\det(A_T,X,Y(X)-Y_S)\ge0. \tag{S92c}
\]

Here is a direct finite-dimensional proof. $L_S$ is linear between
the thresholds $a,b,c$. At $S\le a$ it vanishes. At $S=b$ it is
$(b-a)(bY(c)-cY(b))\ge0$, since $Y(x)/x$ is nondecreasing;
at $S=c$ it is $c\det(\mathbf1,X,Y(X))\ge0$. For $S\ge c$
it is $S\det(\mathbf1,X,Y(X))$. This proves its sign everywhere.

For the other assertions put $H_R(X)=(X-R)_+$,
$D_R=\det(\mathbf1,X,H_R(X))$ and
$E(A,B)=\det(H_A(X),X,H_B(X))$. The nonzero parts of $D_R$ are
$(c-b)(R-a)$ for $a\le R\le b$ and
$(b-a)(c-R)$ for $b\le R\le c$. For $0\le A\le B$,

\[
 E(A,B)=\begin{cases}
 -A D_B,& A\le a,\\
 -a(c-b)(B-A),&a\le A\le B\le b,\\
 -a(b-A)(c-B),&a\le A\le b\le B\le c,\\
 0,& b\le A\ \text{or}\ c\le B.
 \end{cases}
\]

These formulas agree at the shared boundaries and cover all cases.
They show $D_{S,T}=E(S,T)\le0$. By linearity in $Y$, it remains
to prove (S92c) for a linear function and for $Y=H_R$, $R\ge0$.
The linear part cancels. For the hinge the expression is

\[
 (S+T)D_R+E(S,\max(T,R))+E(T,\max(S,R)).
\]

If $R\le S$ this is $(S+T)D_R\ge0$. If $T\le R$ it is
$L_S(H_R)+L_T(H_R)\ge0$. For $S\le R\le T$, set
$B=E(S,T)-E(R,T)$. When $D_R=0$ the displayed formulas give
$B\ge0$. Otherwise they give $B\ge-aD_R$: for
$a\le R\le b$, $R\le T\le b$ this is immediate; for
$a\le R\le b\le T\le c$ first use
$B\ge-a(R-a)(c-T)$; for $b\le R\le T\le c$ use
$B\ge-a(b-a)(c-T)$. In the last two cases $T\ge R$ gives the
claimed bound, and for $T\ge c$ one has $B=0$.
Thus the expression is at least $(S+T-a)D_R\ge0$, since
$D_R>0$ implies $T\ge R>a$. Every convex piecewise linear $Y$
is a linear function plus a nonnegative sum of hinges, completing
the proof. The exact checker also gives an exhaustive coefficient
certificate: all 60 orderings of $a,b,c,S,T,R$ with $a<b<c$ and
$S<T$ have nonnegative cubic gap coefficients, with 240 nonzero
coefficients in total. Equalities follow by continuity.

#### Applying clipping to the polynomial residual

Use the Helmert coordinates divided by their positive factors
$\alpha_j$. Write

\[
 X_j=\sum_{i=1}^j m_i-jm_{j+1},\qquad
 Y_j=-\sum_{i=1}^j m_i^2+jm_{j+1}^2,\qquad m_1=m_2=0.
\]

Starting with $(X_1,Y_1)=(0,0)$, these vertices form an increasing
convex piecewise linear curve $Y(X)$: the first slope is $-m_3>0$,
and the later slopes are $-(m_{j+1}+m_{j+2})$, strictly increasing.
Extend its final segment linearly if necessary. For a hinge
$h_l(t)=(t-m_l)_+$ put $S_l=X_{l-1}$. Its $F$ coordinates are
$A_l=\min(X,S_l)$, its scalar sum is $S_l$, and its $T$ coordinates
are $-Y_{S_l}-m_l A_l$. These identities follow by separating rows
whose endpoint is above or below $m_l$.

For a fixed triple of Helmert indices, (S92b) therefore has linear
hinge term $L_l=\det(V_{S_l},X,Y)$, already nonnegative.
Let $Q_{lj}$ be its symmetric quadratic polarization on two hinges,
and set $C_{lj}=(L_l+L_j)/2+Q_{lj}$. For $l\le j$, direct expansion
gives

\[
 2C_{lj}=C_0+(m_j-m_l)D_{S_l,S_j}\ge0,
\]

using (S92c), $m_j\le m_l$ and $D_{S_l,S_j}\le0$.
The omitted factors $\alpha_{j_1}\alpha_{j_2}\alpha_{j_3}$ are positive.

On the finite set of means write the convex function exactly as
$f=c+\sum_{l=3}^k\eta_lh_l$, where $c=f(m_k)>0$,
$\eta_l\ge0$ and $s=\sum_l\eta_l\le1$. Indeed the $\eta_l$
are successive differences of the ordered secant slopes, with a
final zero slope appended; their sum is the top secant slope.
The constant has $F=0$ and $T=cz$, so all its quadratic cross terms
vanish. Thus each residual coordinate is

\[
 R(f)=kc\,W_3+(1-s)\sum_l\eta_lL_l
                         +\sum_{l,j}\eta_l\eta_jC_{lj}\ge0.
 \tag{S92d}
\]

This proves the comparison for every such convex function, without
an upper bound on the number of means.

#### Mills convexity and the full operator bound

The needed Mills facts are classical; the primary
[Barrow--Cohen paper](https://www.jstor.org/stable/2236746),
*Ann. Math. Statist.* 25 (1954), 405--408, states the equivalent
truncated-moment inequality in (2). For completeness they follow
directly here. If $Z$ is standard normal and
$T=v-Z\mid Z\le v$, let $r=\mathbb ET=v+\lambda(v)>0$ and
$V=\operatorname{Var}T=1-\lambda(v)r>0$.
Differentiation gives $\lambda'=-\lambda r\in(-1,0)$ and
$\lambda''=\lambda(r^2-V)$. The survival function of $T$ is
$S_T(t)=\Phi(v-t)/\Phi(v)$ for $t\ge0$, with strictly increasing
hazard $\lambda(v-t)$. Therefore
$S_T(s+t)<S_T(s)S_T(t)$ for $s,t>0$. Integration yields
$\mathbb ET^2<2(\mathbb ET)^2$, so $V<r^2$ and $\lambda''>0$.
All integrals converge by Gaussian tails. Hence $f(t)=\lambda(u-t)$
is positive, increasing and convex, with derivative in $(0,1)$,
and (S92d) applies at every truncation level.

Multiply the conditional comparison by its positive common density
and integrate over $u$. This proves (S91f) at $p=3$, and its exact
scaling proves (S91e) at every correlation. The positive field and
nonnegative folded kernel from (S91) then give
$\|a_\rho^{-1/2}\mathsf K_3a_\rho^{-1/2}\|\le\rho^2$ on the
entire multiplicity space. The local bound
$a_\rho/\rho\ge\ell/\sigma$ proves the active-space gap
$\ell(1-\rho^2)/\sigma=\ell\sigma$. Unused Hermite degree $r$
multiplies the interaction by $\rho^r$, preserving this bound in
every $d\ge k-1$. This completes (S92).

The exact replay checks the complete clipped-curve polynomial
certificate, finite rational hinge/residual identities and conditional
wedge controls. The all-size hinge decomposition, Mills argument,
Gaussian integrations and complete Hilbert-space spectral passage
remain written mathematics with internal review only.

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
  Every nonstandard two-row block has the strict gap (S87), extending
  the $(k-2,2)$ case (S60), and every sign block is
  settled by (S61)--(S66). The third and fourth exterior hooks are
  settled by (S73) and (S92), respectively.
  The [four-cell](TETRAHEDRAL_STANDARD_GAP.md) and
  [five-cell](FIVE_CELL_A3_BRIDGE.md) comparisons settle (S51) in
  those two cases and complete their full interface gaps. For $k\ge6$
  other blocks with at least three rows also need estimates outside
  the explicit Hermite-floor ranges (S80)--(S81). At fixed correlation
  those ranges cover every sufficiently deep representation, leaving
  only partitions with $O_\rho(\log k)$ boxes outside their first row.
  Equations (S88)--(S90) give the exact antisymmetry hierarchy and a
  common scalar spectral reduction for every fixed representation height.
  The highest projection is bounded below on equivariant fields, and
  (S90) specifies its closed physical image and reconstruction.
  For the remaining exterior hooks, (S91) gives a positive harmonic
  comparison field and exact one-dimensional Gaussian moment formulas.
  Its degree-$3$ comparison is settled by (S92); the sufficient
  intermediate-degree inequalities (S91e) for $p\ge4$ remain open.
  The required scalar spectral-radius and fixed-space bounds are open;
  the next simple odd-kernel positivity argument is explicitly false.
  The exact tied-height reduction (S84) supplies a two-variable
  diagnostic kernel with its original mass constraint; its restricted
  sign problem is also open and would not by itself settle (S51).
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
immanant matrices for every partition through $k=5$ and eight partitions
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
For (S80)--(S83) it independently computes graded character multiplicities
through nine cells and checks rational relative-error identities, the
first two-column instance and rejection of a shallow block outside
the estimate. The unrestricted branching and analytic gap proof are
written arguments; the checker does not certify their hypotheses or
an interface realization of the first possible polynomial degree.
For (S82) it checks the two-row first-derivative transport on a basis
through twelve cells. The conditional covariance theorem, multiplicity
argument and uniform initial-interval deduction remain written proofs.
For the explicit endpoint in (S83), bounded Hermite-power inequalities
and rational endpoint substitutions are checked, with controls that
reject inclusion of an unbalanced constant term or an unsubtracted
first-moment term. The unrestricted nonnegative series comparison
and operator-norm conclusion are the written argument.
For (S84), exact conditional Gaussian matrices at four rational
correlations and formal surface/differentiation coefficient controls
through fourteen cells reject independent adjacent normals and an
omitted adjacent term. They do not certify the surface integrals,
compact compression, trace decomposition or any standard-sector sign.
For (S85)--(S87), the finite Young models check the triple-sign
projection and scalar adjacent corner for two-row shapes. Exact
exponential comparisons on selected half-facet points and the
exceptional-pair factorization test the one-reflection argument
without imposing a full order on the remaining labels. A negative
two-point adjacent-kernel determinant rejects the termwise operator
ordering shortcut. These checks do not certify the all-$k$ branching,
compact supersolution or nontrivial-type spectral transfer proofs.
For (S88)--(S90), formal coset/corner identities are checked in finite
groups through $S_6$, frame/content identities in the finite Young
models, and an exact negative five-cell scalar-reflection witness.
The same models check conjugate highest-sign frame determinants and
commutation, with vanishing above the representation height. Every
vertical two-strip through ten cells is extended to each permitted
size, with a compatible reverse-tableau chain. The all-size coset and
strict frame proofs, imported Young/Pieri/content inputs, reconstruction,
closed-image intertwining and full spectral inclusion are written arguments.
For (S91), the exact controls check positive polynomial fields, their
transported Vandermonde expansion and free Gaussian moments in every
exterior degree for four through eight cells. Finite product measures
check the determinant of unnormalized truncated moments, including
spectator factors, and a formal identity checks the Gaussian moment
recurrence through order nine. These controls do not prove the
remaining intermediate-degree truncated supersolution comparisons (S91e).
For (S92), all 60 real-variable order types for the clipped-curve
inequality are expanded as cubic polynomials in nonnegative successive
gaps, with all 240 nonzero coefficients positive. Ten analogous types
verify the clipped-pair sign. Finite rational hinge and conditional
product controls test the residual identities, including rejection of
an omitted tent term, reversed thresholds and a missing slope bound.
The unrestricted clipping and convex-slope arguments, self-contained
Mills proof, Gaussian integration and full-block operator conclusion
are written mathematics; no external review or formal verification
is inferred from the exact polynomial replay.
The analytic operator reduction, OU parity spectra, tensor decomposition,
dimension lifting, determinant positivity and all-$k$ representation
argument still require written proof review. The optional
quadrature uses NumPy 1.26.4 to check the displayed gap constants and
extra-degree-one block at four chosen correlations, with two grids on
$[0,10]$, and checks the separate perimeter expansion at three angles.
It has no certified tail or discretization errors and proves no
inequality. Its agreement is separate from the exact arithmetic.
