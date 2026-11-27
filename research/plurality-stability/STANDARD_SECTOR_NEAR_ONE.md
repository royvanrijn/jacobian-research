# The Gaussian simplex standard sector near correlation one

At each fixed cell count \(k\ge4\), the ordinary-balanced standard
sector of the simplex interface form has a positive gap for all
correlations sufficiently close to one. Its limiting gap is governed
by the Gaussian Neumann problem on a single facet. Compactness on the
unbounded facet is essential here. The conjugated same-facet Mehler
kernel supplies a Gaussian confinement bound, so a possible sequence
of low-energy modes cannot escape along the facet or concentrate at its
edges. The normalized limiting gap is exactly one at every cell
count. The same compactness argument gives convergence of every
fixed variational level.

## Statement and normalization

Let \(E=\mathbf 1^\perp\subset\mathbb R^k\) and let \(H\subset E\)
be the plane \(x_1=x_2\). Write \(n=k-2=\dim H\), let \(\gamma_H\)
be its standard Gaussian probability measure, and put
\[
 F=\{x\in H:x_2\ge x_l,\ 3\le l\le k\},\qquad
 d\nu=\phi_0\,d\gamma_H|_F,\qquad \phi_0=(2\pi)^{-1/2}.
\]
The standard-sector multiplicity space consists of functions \(h\)
on \(F\) invariant under permutations of labels \(3,\ldots,k\).
The balance constraint is \(\int_F h\,d\nu=0\).
For such \(h\), use the scalar form from
[(S49)–(S51)](SIMPLEX_INTERFACE_SYMMETRY.md#every-cell-count-the-standard-sector-is-a-scalar-positive-kernel):
\[
 q_\rho^{\rm std}(h)=\frac1\rho\left[
   \int_F a_\rho h^2\,d\nu-
   \iint_{F^2}K_{{\rm std},\rho}(x,y)h(x)h(y)\,d\nu(x)d\nu(y)
 \right].
\]
In the tensor realization \(f^{v,h}\) of (S49), the full interface
form and the facet norm squared are respectively \(k|v|^2 q_\rho^{\rm std}(h)\)
and \(k|v|^2\|h\|_{L^2(\nu)}^2\). Consequently the scalar quotient is
also the standard-sector gap in the full facet norm.

Put \(\sigma_\rho=\sqrt{1-\rho^2}\) and denote the first invariant
Gaussian Neumann Rayleigh value by
\[
 \lambda_k=\inf_{\substack{0\ne h\in W^{1,2}(F^\circ,\gamma_H)\\
             h\ {\rm invariant},\ \int_Fh\,d\gamma_H=0}}
                  \frac{\int_F|\nabla h|^2\,d\gamma_H}
                       {\int_Fh^2\,d\gamma_H}.
\]
Here \(F^\circ\) means interior in \(H\), and the boundary has
\(\gamma_H\)-measure zero. Denote the balanced scalar interface gap by
\(g_k(\rho)=\inf q_\rho^{\rm std}(h)/\|h\|_{L^2(\nu)}^2\), with the
infimum over nonzero invariant balanced \(L^2(\nu)\) functions.

**Theorem.** For each fixed \(k\ge4\),
\[
 \lim_{\rho\uparrow1}\frac{\rho\,\sigma_\rho}
                     {2\phi_0(1-\rho)}\,g_k(\rho)
              =\lambda_k-1>0 .
 \tag{SN1}
\]
For every \(k\ge4\), \(\lambda_k=2\), with eigenspace spanned by
\(|x|^2-n\). Thus the constant on the right is exactly one.
For \(k=5\), the standard result and
[(S76)](SIMPLEX_INTERFACE_SYMMETRY.md#strict-supersolution-and-the-five-cell-reduction)
give the full balanced interface gap off rotations near one, and
[(S42)](SIMPLEX_INTERFACE_SYMMETRY.md#every-cell-count-the-active-gap-lifts-to-all-ambient-dimensions)
gives a gap independent of the ambient dimension \(d\ge4\) at each
fixed correlation in that interval. The cutoff in this argument
depends on \(k\) and is qualitative. The remaining symmetry sectors
for \(k\ge6\) require separate estimates. This is an interface
theorem on the simplex, and comparison with distant balanced
partitions is a separate question.

## Same-facet reduction and the two-sided test

Let \(m_\rho(x,y)\) be the Mehler kernel of \(\gamma_H\) and put
\(w_l=(e_1+e_2)/2-e_l\) for \(3\le l\le k\). On \(H\), the compression
of \((2l)\) is \(I-w_lw_l^T\), so the planar kernels of the
adjacent-facet terms are
\[
 m_{\rho,l}(x,y)=m_\rho(x,y)
              e^{-\rho(w_l\cdot x)(w_l\cdot y)/\sigma_\rho^2}.
\]
All these kernels are nonnegative, and \(m_{\rho,l}\le m_\rho\) on
\(F^2\). The expression (S49) in the ambient Gaussian has
\(K_{{\rm std},\rho}=\sigma_\rho^{-1}(2m_\rho+\sum_l m_{\rho,l})\)
since the ambient Mehler prefactor has one extra factor
\(\sigma_\rho^{-1}\), and the transpositions preserve the squared
norms in the exponent. For \(m=m_\rho\) or \(m_{\rho,l}\), write
\[
 J_m(h)=\frac12\iint_{F^2}m(x,y)(h(x)-h(y))^2\,d\gamma_H(x)d\gamma_H(y).
\]
Use \(J_\rho=J_{m_\rho}\), \(J_{\rho,l}=J_{m_{\rho,l}}\), and
\(b_\rho(x)=\sigma_\rho a_\rho(x)/(2\phi_0\rho)\). The row identity
\(\int_F K_{{\rm std},\rho}(x,y)\,d\nu(y)=a_\rho(x)/\rho\) then gives
\[
 \frac{\rho\sigma_\rho}{2\phi_0(1-\rho)}
       \frac{q_\rho^{\rm std}(h)}{\|h\|_{L^2(\nu)}^2}
 =\frac{J_\rho(h)+\tfrac12\sum_lJ_{\rho,l}(h)}
       {(1-\rho)\|h\|_{L^2(F,\gamma_H)}^2}
   -\frac{\int_F b_\rho h^2\,d\gamma_H}{\|h\|_{L^2(F,\gamma_H)}^2}.
 \tag{SN2}
\]
In [(B7)](POSITIVE_CORRELATION_CONTINUATION.md#the-full-interface-form-at-a-fixed-positive-correlation),
\(0<b_\rho\le1\). Its cone formula also gives
\(b_\rho(x)\to1\) for every \(x\in F^\circ\). Indeed the translated
region \(T_1-(\rho/\sigma_\rho)x\) in that formula converges pointwise
off boundaries to the halfspace \(z_1\ge z_2\). Each other inequality
becomes automatic for fixed \(z\) because \(x_1>x_l\), \(l\ge3\).
Integrating the positive normal coordinate against the full Gaussian
on that halfspace gives \(\phi_0\). This proves the limit by dominated
convergence.

For a fixed smooth \(h\) on \(H\) with bounded gradient, substitute
\(y=\rho x+\sigma_\rho z\) in \(J_\rho\) and \(J_{\rho,l}\). The
indicator \(1_F(\rho x+\sigma_\rho z)\) tends to one for \(x\in F^\circ\)
and each fixed \(z\). The quotient
\((h(\rho x+\sigma_\rho z)-h(x))/\sigma_\rho\) tends to
\(z\cdot\nabla h(x)\), and its square is bounded by a constant times
\((|z|+|x|)^2\). In the adjacent terms the additional exponential
tends to zero for such \(x,z\). Since
\(\sigma_\rho^2/(2(1-\rho))\to1\), dominated convergence proves
\[
 \frac{J_\rho(h)}{1-\rho}\longrightarrow\int_F|\nabla h|^2d\gamma_H,
 \qquad
 \frac{J_{\rho,l}(h)}{1-\rho}\longrightarrow0 .
 \tag{SN3}
\]
These tests are dense in the invariant balanced Gaussian Sobolev space.
For an explicit approximation, first cut off radially in
\(W^{1,2}(F,\gamma_H)\). For a compactly supported cutoff \(g\), choose
\(b\in F^\circ\) and \(c>0\) so \(b+B_c\subset F^\circ\).
Mollify the zero extension of \(g\) at \(x+\epsilon b\) at scale
\(c\epsilon/2\). For every \(x\in F\) the mollification uses points of
\(F^\circ\). Its gradient there is the same mollification of the
zero extension of the weak gradient of \(g\), by the local
integration by parts identity. Strong \(L^2(H)\) continuity of
translations and mollification for these compactly supported data
gives strong \(W^{1,2}(F,\gamma_H)\) approximation.
Finally average over the facet stabilizer and subtract the ordinary
Gaussian mean on \(F\). These operations preserve convergence, and
the resulting tests are smooth with bounded gradient. Thus (SN2)
and (SN3) give the upper bound \(\lambda_k-1\) for the limiting
rescaled infimum.

## Gaussian confinement and compactness

The sharp lower bound requires estimates for functions varying with
\(\rho\). Put \(t=1-\rho\) and let \(\phi_H\) denote the Gaussian
density. Conjugate by \(u=\sqrt{\phi_H}\,h\) and introduce the
Lebesgue kernel \(\widetilde m_\rho(x,y)=
\sqrt{\phi_H(x)\phi_H(y)}\,m_\rho(x,y)\).
Let \(m^F_\rho(x)=\int_Fm_\rho(x,y)\,d\gamma_H(y)\) and
\(\widetilde m^F_\rho(x)=\int_F\widetilde m_\rho(x,y)\,dy\).
Expanding the two squared differences gives
\[
 J_\rho(h)=\frac12\iint_{F^2}\widetilde m_\rho(x,y)
                       (u(x)-u(y))^2\,dx\,dy
          +\int_F(m_\rho^F-\widetilde m_\rho^F)u^2\,dx .
 \tag{SN4}
\]
The exact exponent of the conjugated kernel is
\[
 \widetilde m_\rho(x,y)
 =\frac{1}{(2\pi\sigma_\rho^2)^{n/2}}
   \exp\left[-\frac{1+\rho}{8(1-\rho)}|x-y|^2
             -\frac{1-\rho}{8(1+\rho)}|x+y|^2\right].
\]
For \(Z\sim\gamma_H\), direct Gaussian integration and completion of
the square give
\[
 \begin{split}
 m_\rho^F(x)&=\Pr[Z+(\rho/\sigma_\rho)x\in F],\\
 \widetilde m_\rho^F(x)
 &=s_\rho^n e^{-\sigma_\rho^2|x|^2/(4(1+\rho^2))}
       \Pr[Z+s_\rho(\rho/\sigma_\rho)x\in F],
 \qquad s_\rho=\sqrt{2/(1+\rho^2)} .
 \end{split}
\]
Here \(m_\rho^F\ge p_F:=\gamma_H(F)>0\) by the cone property and
\(x\in F\). The relative threshold increase from \(1\) to \(s_\rho\)
raises the probability by at most \(C_k t\), uniformly on \(F\).
For each face normal directed inward and normalized to unit length,
the original threshold is \(-a\) with \(a\ge0\) and the moved one
is \(-s_\rho a\). A union bound over the face intervals bounds the
difference of the two cone probabilities by the sum of
\(\Phi(s_\rho a)-\Phi(a)\le(s_\rho-1)\sup_{v\ge0}v\phi_1(v)\).
Since \(s_\rho^n=1+O_k(t)\) and
\(\sigma_\rho^2/(4(1+\rho^2))\ge t/4\), we obtain
\[
 m_\rho^F(x)-\widetilde m_\rho^F(x)
    \ \ge\ p_F(1-e^{-t|x|^2/4})-C_k t .
 \tag{SN5}
\]
This global bound controls the Gaussian \(L^2\) tails and prevents
possible concentration at infinity from being lost in a bounded-domain
limit.

For a sequence \(\rho_j\uparrow1\) and \(h_j\) with
\(\|h_j\|_{L^2(F,\gamma_H)}=1\) and \(J_{\rho_j}(h_j)\le C(1-\rho_j)\),
(SN4)–(SN5) show that the conjugated difference integral is
\(O_k(1-\rho_j)\) and
\[
 \limsup_j\int_{F\setminus B_R}|u_j|^2\,dx\ \le\ C'/R^2
 \quad(R>0).
 \tag{SN6}
\]
For the latter, bound \(1-e^{-t|x|^2/4}\) below by \(1-e^{-tR^2/4}\)
outside \(B_R\) and take the limit of \(t/(1-e^{-tR^2/4})\).

On the bounded domain \(\Omega_R=F^\circ\cap B_R\),
\(\widetilde m_{\rho_j}(x,y)\ge c_R\sigma_{\rho_j}^{-n}\) whenever
\(|x-y|\le\sigma_{\rho_j}\). We use
[Augusto C. Ponce, *An estimate in the spirit of Poincaré's inequality*,
J. Eur. Math. Soc. 6 (2004), Theorem 1.2, p. 2](sources/ponce-nonlocal-poincare-2004.pdf).
The theorem gives \(L^2\) relative compactness for bounded
\(L^2\) sequences on bounded Lipschitz domains in dimension \(n\ge2\)
with bounded difference integrals whose kernel is
\(\eta_j(x-y)/|x-y|^2\), where \(\eta_j\) are radial probability
densities concentrating at zero. Here \(\Omega_R\) is bounded,
convex, with nonempty interior, and hence is Lipschitz.
Take
\(\eta_j(z)=c_n\sigma_{\rho_j}^{-n-2}|z|^2 1_{\{|z|\le\sigma_{\rho_j}\}}\),
with \(c_n\) normalizing its integral. The small-distance
kernel bound and the \(O_k(1-\rho_j)\) difference integral give
the required bound because \(1-\rho_j=O(\sigma_{\rho_j}^2)\).
A subsequence converges strongly in \(L^2(\Omega_R)\) for each \(R\)
by a diagonal argument. The tails (SN6) give strong convergence in
\(L^2(F,dx)\). By conjugation, \(h_j\to h\) strongly in
\(L^2(F,\gamma_H)\). The limit has unit norm and preserves balance
and the facet symmetry.

## Sharp lower limit

For \(U\Subset F^\circ\) and \(L>0\), use instead
\(y=x+\sigma_{\rho_j}z\) in the unconjugated same-facet term.
For all sufficiently large \(j\), \(U+\sigma_{\rho_j}B_L\subset F^\circ\).
Consequently \(J_{\rho_j}(h_j)/(1-\rho_j)\) is bounded below by
\[
 \frac{\sigma_{\rho_j}^2}{2(1-\rho_j)}
     \int_U\int_{B_L}
       \left|\frac{h_j(x+\sigma_{\rho_j}z)-h_j(x)}
                        {\sigma_{\rho_j}}\right|^2
       \phi_H\!\left(z+\frac{1-\rho_j}{\sigma_{\rho_j}}x\right)
      \,dz\,d\gamma_H(x).
 \tag{SN7}
\]
The weights in this compact region are bounded below and converge
uniformly to the unshifted Gaussian weights. The difference
quotients therefore have weak \(L^2(U\times B_L)\) subsequential
limits. If \(G\) is such a limit, change variables in the translated
term of the pairing with a smooth compactly supported test \(\psi\).
The strong local \(L^2\) convergence of \(h_j\) gives
\[
 \int_{U\times B_L}G(x,z)\psi(x,z)\,dx\,dz
     =-\int_{U\times B_L} h(x) z\cdot\nabla_x\psi(x,z)\,dx\,dz .
\]
Testing with \(\psi(x,z)=\chi(x)z_i\eta(z)\), where \(\eta\) is a
radial test with \(\int z_i^2\eta(z)\,dz=1\), identifies
\(\partial_i h\in L^2(U)\). The same identity for general tests
then gives \(G=z\cdot\nabla h(x)\).
Weak lower semicontinuity in (SN7), followed by exhaustion in \(U\)
and \(L\), proves
\[
 \liminf_j\frac{J_{\rho_j}(h_j)}{1-\rho_j}
       \ \ge\ \int_F|\nabla h|^2\,d\gamma_H .
 \tag{SN8}
\]
The constant follows from the limiting factor
\(\sigma_{\rho_j}^2/(2(1-\rho_j))\to1\) and the identity Gaussian
covariance in \(z\). In particular the limiting gradient is globally
square integrable.

Apply this to approximate unit-norm minimizers along a sequence on
which the rescaled gap approaches its lower limit. The fixed tests
(SN3) bound these rescaled quotients above and (SN2) bounds them below
by \(J_{\rho_j}(h_j)/(1-\rho_j)-1\). This supplies the energy
bound used in the compactness proof. By balance and symmetry of the
unit-norm limit, (SN8) and (SN2) prove a lower bound of
\(\lambda_k-1\). This is the reverse inequality to the fixed-test
upper bound.

## The invariant Gaussian Neumann gap

Split \(h\) in polar coordinates \((r,\theta)\) on \(F\) into its
angular mean \(h_0(r)\) and the angular-mean-zero remainder
\(h_1(r,\theta)\). The two terms are orthogonal in norm and energy
for the Gaussian measure, whose radial density is proportional to
\(w(r)=r^{n-1}e^{-r^2/2}\). The Gaussian radial Rayleigh value on
balanced functions is \(2\). Indeed the radial extension to the
entire \(H\) is orthogonal to constants and all linear functions,
so full-space Hermite expansion gives the bound. For radial
Sobolev functions the extension across the origin can be justified
by making \(h_0\) constant below a sequence \(r_i\downarrow0\)
with \(r_i^n|h_0(r_i)|^2\to0\). Such a sequence follows from
integrability in the radial norm, and the modifications converge
in the radial norm and energy. The function \(r^2-n\) attains \(2\).

For \(n\ge3\), let \(S=F^\circ\cap S^{n-1}\) and take its first
nonconstant Neumann eigenvalue \(\Lambda\) for spherical area.
We use [Andrews–Clutterbuck–Hauer, *Non-concavity of the Robin ground
state*, Cambridge J. Math. 8 (2020), Theorem 9.1, pp. 296–297](sources/andrews-clutterbuck-hauer-robin-ground-state-2020.pdf).
The theorem applies to convex open spherical domains in \(S^{n-1}\)
for \(n\ge3\), including nonsmooth domains. It gives
\(\Lambda\ge n-1\), with equality exactly when the cone over the
domain has a linear factor. Here \(S\) is a geodesically convex
intersection of hemispheres with nonempty interior. The normals
\(w_l\) have Gram matrix \(I+\tfrac12\mathbf1\mathbf1^T\) and span
\(H\), so the cone has trivial lineality. Thus \(\Lambda>n-1\).
For \(n=2\) and \(k=4\), invariance identifies the angular problem
with the interval \((0,\alpha)\), \(\alpha=\arctan\sqrt2\), by
reflection in the bisector of the full facet. The Neumann value
in that interval is \(\Lambda=(\pi/\alpha)^2\) by cosine expansion.

Let \(\beta>0\) solve \(\beta(\beta+n-2)=\Lambda\), where
\(\Lambda\) can be any lower bound for the angular Neumann gap
on the invariant space. The angular
Poincaré inequality lower bounds the energy of \(h_1\) by the
radial-derivative energy plus \(\Lambda\int_0^\infty
\|h_1(r,\cdot)\|^2 r^{-2}w(r)\,dr\). In the angular \(L^2\) norm,
for \(g=h_1\), integration by parts gives
\[
 \int_0^\infty
   \big(\|g'\|^2+\Lambda r^{-2}\|g\|^2\big)w(r)\,dr
   =\beta\int_0^\infty\|g\|^2w(r)\,dr
     +\int_0^\infty\|g'-\beta g/r\|^2w(r)\,dr .
\]
For compact radial support away from zero this follows from
\((w/r)'=((n-2)/r^2-1)w\). To pass to all finite-energy \(g\), the
angular Poincaré inequality supplies \(g/r\in L^2(w\,dr)\).
Use radial cutoffs near zero and infinity with derivative bounded
by a constant times \(r^{-1}\) on their transition annuli. The
resulting gradient errors vanish by integrability of \(g/r\).

For \(n\ge3\), the estimate \(\Lambda>n-1\) already gives
\(\beta>1\) for unrestricted facet functions. To identify the
standard-sector constant, use the sharper estimate on the ordered
facet.

### The ordered facet has angular Neumann gap strictly above \(2n\)

In the coordinates \(g_i=x_2-x_{i+2}\), \(1\le i\le n\), order the
cone by \(0\le g_1\le\cdots\le g_n\). Restriction identifies the
invariant Sobolev space with the ordinary Sobolev space on this
chamber. Conversely, extend by permuting the untied labels. Traces
match across the ordering walls since each transposition fixes
its wall pointwise. The full-facet squared norm and Dirichlet
energy are both \(n!\) times their chamber values. The same
argument applies on the spherical link.

In these coordinates the Euclidean metric is
\((I+J/2)^{-1}\), with \(J=\mathbf 1\mathbf 1^\top\), since the
gradients of the \(g_i\) have Gram matrix \(I+J/2\). The inward
wall covectors are
\[
 b_0=e_1,\qquad b_j=e_{j+1}-e_j\quad(1\le j<n).
\]
They are independent, and their pairwise dual-metric inner
products are \(-1\) for consecutive indices and \(0\) for all
other distinct pairs. Indeed each \(b_j\), \(j\ge1\), has sum
zero, so the \(J/2\) term contributes only to \(\|b_0\|^2\).
Thus the ordered spherical chamber is a simplex with all
interior dihedral angles at most \(\pi/2\).

For \(n\ge3\), Corollary 1.3 of Edelen–Li gives a first
nonconstant Neumann eigenvalue at least \(2n\) on this chamber
in \(S^{n-1}\). In the equality case an eigenfunction
is the restriction of a homogeneous harmonic quadratic.
[Edelen–Li, *Sharp Neumann eigenvalue estimates and \(C^2\)
elliptic regularity in non-obtuse polyhedral domains*, arXiv
v1, Corollary 1.3](sources/edelen-li-2608.17194v1.pdf)

Suppose that quadratic is \(x^\top Mx\) with \(M=M^\top\)
and \(\operatorname{tr}M=0\). Neumann data on the relative
interior of each face give \(x^\top M\nu=0\) there, since
the unit normal \(\nu\) to the wall already lies in the
tangent space of the sphere along the wall. The face spans
the wall plane, so \(M\nu\) is proportional to \(\nu\).
Consecutive wall normals have nonzero inner product.
Eigenvectors of a real symmetric matrix for different
eigenvalues are orthogonal, so all wall normals have the
same eigenvalue. They form a basis, so \(M\) is scalar.
Harmonicity gives \(M=0\), excluding equality. For \(n=2\),
the explicit interval value above gives the same strict
inequality. In particular,
\[
 \Lambda_{\mathrm{ord}}>2n,\qquad
 \lambda_k=2 \qquad(k\ge4).
 \tag{SN10}
\]
In the radial separation this strict inequality gives
\(\beta>2\); equality in the Gaussian Rayleigh quotient
therefore requires a radial function. Full-space Hermite
expansion leaves precisely the span of \(r^2-n\).

The same argument covers chains of tied blocks. Let \(B_1,\ldots,B_s\)
be disjoint coordinate blocks of positive sizes \(m_1,\ldots,m_s\)
with \(s\ge2\). Restrict to coordinates that take the constant
value \(t_i\) on \(B_i\), satisfy \(\sum_i m_i t_i=0\), and are
ordered \(t_1\ge\cdots\ge t_s\). The \(s-1\) wall normals are
\[
 n_i=\frac{\mathbf1_{B_i}}{m_i}
       -\frac{\mathbf1_{B_{i+1}}}{m_{i+1}},\qquad 1\le i<s.
\]
They form a basis of the centered block-constant subspace.
Consecutive inner products are \(-1/m_{i+1}\), and all other
distinct pairs are orthogonal. Thus for \(s\ge4\) the same
Neumann comparison and connectedness argument give a strict
angular gap above \(2(s-1)\). For \(s=3\) there is an acute
interval with gap above \(4\). For \(s=2\) the chamber is a
half-line, and even extension gives the same Gaussian
Rayleigh value \(2\). In every case the first ordinary Gaussian
Neumann eigenspace is exactly the span of \(|x|^2-(s-1)\).

### Explicit angular bounds in five and six cells

For the invariant functions we can use the ordered sector of \(F\)
where \(0\le g_1\le\cdots\le g_n\), with
\(g_i=x_2-x_{i+2}\). In the reference space
\(E'=\mathbf1^\perp\subset\mathbb R^{n+1}\) with coordinates
\(y_0,\ldots,y_n\), the chamber \(y_0\ge y_1\ge\cdots\ge y_n\)
has the same parameters \(g_i=y_0-y_i\). The squared norms in the
two spaces are
\[
 \|x\|^2=\sum_i g_i^2-\frac{(\sum_i g_i)^2}{n+2},\qquad
 \|y\|^2=\sum_i g_i^2-\frac{(\sum_i g_i)^2}{n+1}.
\]
The linear map from \(y\) to \(x\) therefore preserves lengths
orthogonal to the all-equal-gap axis and dilates that axis by
\(a=\sqrt{2(n+1)/(n+2)}\).

Let \(C\) and \(A\) be the spherical chambers of the reference and
actual cones. If \(y\in C\) has unit norm, its coordinate along
the positive axis is \(X=(\sum_i g_i)/\sqrt{n(n+1)}\). Since
\(\sum_i g_i^2\le(\sum_i g_i)^2\), one has \(X\ge1/n\).
Set
\[
 s=\sqrt{1+(a^2-1)X^2},\qquad s_*=\sqrt{1+(a^2-1)/n^2}.
\]
The map to \(A\) applies the linear dilation and divides by \(s\).
Writing \(X=\cos\theta\), the new angle \(\theta'\) satisfies
\(\sin\theta'=\sin\theta/s\) and
\(d\theta'/d\theta=a/s^2\). The metric factors on tangent lengths
are therefore \(a/s^2\) in the meridian direction and \(1/s\)
in each transverse direction, with area Jacobian \(a/s^n\).

The Neumann gap of \(C\) is at least \(3(n+1)\).
Extend functions on \(C\) to the sphere by permutation symmetry.
The traces coincide at the chamber walls because the adjacent
transposition fixes its wall pointwise. The reflected pieces thus
give a spherical \(H^1\) function, with both the norm and energy
multiplied by \((n+1)!\). In the spherical harmonic expansion
of a mean-zero invariant function, degrees one and two vanish.
The only invariant linear functional on \(E'\) is zero, and
an invariant quadratic is a multiple of \(|y|^2\), whose harmonic
part is zero. For the quadratic statement one can first extend
to \(\mathbb R^{n+1}\) by orthogonal projection onto \(E'\), where
symmetric quadratics are combinations of \(\sum_i y_i^2\) and
\((\sum_i y_i)^2\). The degree-three spherical eigenvalue on
\(S^{n-1}\) is \(3(n+1)\).

Let \(U\) be the pullback to \(C\) of a function on \(A\).
Because \(s\le a\), the smaller of the two gradient multipliers
squared is \(s^4/a^2\). The energy on \(A\) is thus at least
\(a^{-1}\min_C s^{4-n}\) times the energy of \(U\) on \(C\).
For the centered squared norm on \(A\), use the mean of \(U\)
on \(C\) as a test constant in the variance infimum. Its value
is at most \(a/s_*^n\) times the centered squared norm on \(C\).
For \(n\le4\) this proves the angular lower bound
\[
 \Lambda_A\ge 3(n+1)\frac{s_*^4}{a^2}.
 \tag{SN9}
\]
For \(n=3\) the right side is \(128/15>6\).
For \(n=4\) it is \(625/64>8\). Thus in both cases we may use
an angular bound \(\Lambda>2n\), which gives \(\beta>2\) in the
radial separation. Since \(r^2-n\) is invariant and attains the
radial value, \(\lambda_5=\lambda_6=2\).

## Higher limiting levels

Let \(0=\Lambda_0<\Lambda_1\le\Lambda_2\le\cdots\) be the ordered
spherical chamber Neumann values with multiplicities.
Write \(\beta_j\ge0\) for the solution of
\(\beta_j(\beta_j+n-2)=\Lambda_j\), with \(\beta_0=0\). On the
balanced subspace \(X\subset L^2(F,\gamma_H)\) of invariant
functions, let \(\xi_d(\rho)\) be the \(d\)-th min–max value of
the quadratic form in (SN2), starting at \(d=1\).

**Theorem (Higher limiting levels).** For every fixed \(k\ge4\)
and \(d\ge1\), \(\xi_d(\rho)\) converges as \(\rho\uparrow1\) to
the \(d\)-th member, in nondecreasing order with multiplicity,
of
\[
 \big\{\beta_j+2m-1:\ j,m\ge0,\ (j,m)\ne(0,0)\big\}.
 \tag{SN11}
\]
For \(\rho\) sufficiently close to one at each fixed \(k,d\),
these first \(d\) min–max values are discrete eigenvalues.
The angle-independent family contributes \(1,3,5,\ldots\)
to the limit. Its first member precedes every angular value
by (SN10). The ordering at higher levels is given by (SN11).

To obtain the limit spectrum, diagonalize on the spherical
chamber in angular eigenfunctions \(Y_j\). With \(w(r)\) the
radial weight above, the radial operator in \(Y_j\) has form
\(\int(|u'|^2+\Lambda_j r^{-2}|u|^2)w\,dr\). Substitute
\(u=r^{\beta_j}v\) and use the radial square identity. The
form becomes \(\beta_j\|u\|^2\) plus the derivative form in
the weight \(r^{n+2\beta_j-1}e^{-r^2/2}\,dr\). For \(j>0\)
the identity extends from compactly supported radial tests
by the cutoffs already used above. For \(j=0\) there is
no inverse-square term. In \(s=r^2/2\) the derivative form
is \(2\int s|\partial_s v|^2 s^{\alpha_j}e^{-s}\,ds\),
with \(\alpha_j=n/2+\beta_j-1\).

The radial differential operator in that weight is
\(-2[s\partial_s^2+(\alpha_j+1-s)\partial_s]\). It preserves
each polynomial degree filtration and has leading
coefficient \(2m\) on monic polynomials of degree \(m\).
Integration by parts makes it symmetric, so the associated
Laguerre polynomials, obtained by orthogonalizing monomials,
are its eigenfunctions. Their separated functions belong
to the original form domain and have zero radial flux at
zero against smooth tests, since for \(j>0\) the flux is
\(O(r^{n+\beta_j-2})\) and \(\beta_j>2\), and for \(j=0\)
the radial derivative is \(O(r)\). Cutoffs extend the weak
eigenvalue identity to the form domain. The polynomials are
dense in the Gamma-weighted space. Indeed the Laplace
transform of an \(L^2\) function times the Gamma density is
analytic on \(\operatorname{Re}z<1/2\). Orthogonality to every
monomial makes its Taylor series at zero vanish. Uniqueness
of the Fourier transform then gives the zero function.
Thus \(r^{\beta_j}L_m^{(\alpha_j)}(r^2/2)Y_j\) form a
complete eigenbasis for the Gaussian form. The shift by
\(-1\) and the omission of the constant pair give (SN11).
Since \(\beta_j\to\infty\), only finitely many pairs occur
below any fixed energy.

For the upper limit in (SN11), approximate each member of
the first \(d\) limit eigenspaces, taking only \(d\)
orthonormal vectors, by the balanced smooth tests used in
(SN3). Polarization gives convergence of the form matrices
on the resulting finite-dimensional spaces. Taking the
approximation error to zero proves the min–max upper bound.

For the lower limit, take \(\rho_i\uparrow1\) realizing the
lower limit and \(d\)-dimensional spaces with Rayleigh
suprema within \(o(1)\) of \(\xi_d(\rho_i)\). Choose
orthonormal bases. Each basis vector has bounded
\(J_{\rho_i}/(1-\rho_i)\) by (SN2), so one finite diagonal
extraction gives strong limits for the entire basis.
The limit vectors stay orthonormal and balanced. For each
fixed linear combination, (SN8) gives the gradient lower
limit. The potential term converges to minus the squared
norm by strong convergence and \(0<b_\rho\le1\). This
bounds the limiting Rayleigh supremum on the span of the
limit vectors by the lower limit of \(\xi_d(\rho_i)\), as
required.

For fixed \(0<\rho<1\), the integral operator of the
rescaled form is compact. Its kernel is bounded in
absolute value by a finite multiple of the Mehler kernel,
which is Hilbert–Schmidt by the ambient Hermite expansion
with squared-norm sum \((1-\rho^2)^{-n}\).
The diagonal term of (SN2) is bounded below by
\(p_F/(1-\rho)-1\), using the self-facet row bound in (SN5)
and nonnegativity of the adjacent rows. These properties
persist under compression to \(X\). The bottom of the
essential spectrum is therefore at least that lower bound.
It tends to infinity, whereas every fixed variational
level has a finite limit, proving the discreteness
assertion. In particular the limiting first level is simple
and the limiting second level is \(\min\{3,\beta_1-1\}>1\).
The first eigenvalue of the rescaled form is simple for
\(\rho\) sufficiently close to one, and its normalized
eigenfunction converges up to sign to the normalized
\(r^2-n\). This follows by the first two level limits
and the same compactness and lower-limit argument applied
to first eigenfunctions.

External mathematical review and formal verification remain open.
