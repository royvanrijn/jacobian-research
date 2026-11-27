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
permanent formulas on an ordered chamber. These reductions do not supply
the missing active-space gap for
$k\ge4$ at arbitrary correlation. The remaining higher-cell blocks and
global optimality remain open. External review and formal verification
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

For $k\ge4$, (S4) supplies the precise symmetry blocks but no sign for
their remaining infinite dimensional parts. The following reductions
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
remaining sign-sector assertion is

\[
 \left\|a_\rho^{-1/2}\mathcal K_\rho
                         a_\rho^{-1/2}\right\|_{L^2(W,\nu_{12})}<1.
 \tag{S48}
\]

The operator is compact, self-adjoint and nonnegative as a quadratic
form, by restriction of the Mehler form. A positive kernel alone gives
neither this norm bound nor comparison with $a_\rho$. A positive
supersolution $w$ satisfying
$\mathcal K_\rho w\le\theta_\rho a_\rho w$ with
$\theta_\rho<1$ would suffice by the weighted Schur bound; none is
established here at arbitrary correlation.

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
balanced gap for $k\ge4$, including (S48) and the other blocks in (S4),
is still an open mathematical obligation. Neither a global simplex
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
reflection signs, not this weighted Poincare inequality. The higher
hook and other blocks in (S4) require their own estimates as well.

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

* Prove the active-space higher-cell gap in (S4). The sign block has the
  positive determinant kernel (S46), with the precise remaining norm
  bound (S48); the standard block has the positive kernel (S49)--(S50)
  and the unresolved weighted Poincare inequality (S51).
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
The analytic operator reduction, OU parity spectra, tensor decomposition,
dimension lifting, determinant positivity and all-$k$ representation
argument still require written proof review. The optional
quadrature uses NumPy 1.26.4 to check the displayed gap constants and
extra-degree-one block at four chosen correlations, with two grids on
$[0,10]$, and checks the separate perimeter expansion at three angles.
It has no certified tail or discretization errors and proves no
inequality. Its agreement is separate from the exact arithmetic.
