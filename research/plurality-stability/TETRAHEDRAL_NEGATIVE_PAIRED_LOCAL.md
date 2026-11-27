# Local stability of the negative tetrahedral pair

Let \(T\) be the regular tetrahedral partition of \(\mathbb R^d\), \(d\ge3\), extended cylindrically from its three-dimensional span. For measurable \(\Delta_4\)-valued functions \(F,G\), let
\[
 J_{-r}(F,G)=\sum_i\langle F_i,T_{-r}G_i\rangle_{\gamma_d},
 \qquad
 d_T(F,G)=\inf_{Q\in O(d)}
 \sum_i\big(\|F_i-1_{QT_i}\|_1+\|G_i-1_{QT_i}\|_1\big),
\]
where \(\Delta_4=\{p_i\ge0:\sum_i p_i=1\}\) and \(T_{-r}\) is the Gaussian Ornstein–Uhlenbeck operator. Assume only the corresponding-mass constraints \(\int F_i=\int G_i\) for \(1\le i\le4\). The four common masses are free.

**Theorem.** There exist \(\varepsilon,\eta,c>0\), independent of \(d\ge3\), such that, if \(|r-1/3|<\varepsilon\) and \(d_T(F,G)<\eta\), then
\[
 J_{-r}(F,G)-J_{-r}(1_T,1_T)\ \ge\ c\,d_T(F,G)^2. \tag{N1}
\]
The weighted Hessian at \(r=1/3\) has a positive gap off common rotations, with the same gap bound in every \(d\). Taking \(G=F\) includes arbitrary nearby assignments at negative correlation with free masses. The independent pair problem requires matching corresponding masses. Indeed, opposite translations give a strictly negative Hessian direction when that constraint is removed.

Heilman's [*Hyperstable sets with voting and algorithmic hardness applications*](https://www.stevenheilman.org/~heilman/papers/hyper.pdf), Problem 7.1 and Lemma 7.6, formulate the paired comparison after replacing \(G(x)\) by \(G(-x)\) and give its second variation. Theorem 8.5 there assumes global minimality and hyperstability. The argument below instead evaluates the weighted second variation at the tetrahedron itself, then uses measurable normal flux. The local cost tends to zero along remote facets. The measurable comparison uses exterior estimates in the active coordinates and Gaussian hypercontractivity in the unused coordinates.

## Interface and two exact eigenspaces

Use the active space \(E=\{x\in\mathbb R^4:\sum_i x_i=0\}\). Orient each pair \(e=ij\), \(i<j\), by \(d_e=e_i-e_j\) and \(n_{ij}=d_e/\sqrt2\), pointing into \(T_i\). Let \(\nu_e=\phi_d\,d\mathcal H^{d-1}\) on \(F_{ij}=T_i\cap T_j\), and use signed speeds \(f_{ji}=-f_{ij}\), with positive \(f_{ij}\) denoting outward motion of \(T_i\). Define
\[
 U_{r,ij}=T_{-r}(1_{T_j}-1_{T_i}),\qquad
 a_r|_{F_{ij}}=\partial_{n_{ij}}U_{r,ij},\qquad
 (Cf)_i=\sum_{j\ne i}\int f_{ij}\,d\nu_{ij}.
\]
The kernel acting between these facets is
\[
 K_r(x_e,y_h)=(d_e\cdot d_h) M_{-r}(x,y),\qquad
 M_{-r}(x,y)=\sigma^{-d}e^{-\beta(|x|^2+|y|^2)-\tau x\cdot y},
 \quad
 \sigma^2=1-r^2,\quad \beta=\frac{r^2}{2\sigma^2},\quad \tau=\frac r{\sigma^2}.
\]
The measures \(\lambda_i^f=\sum_{j\ne i}f_{ij}\nu_{ij}\) yield
\(K_r(f,g)=\sum_i\iint M_{-r}\,d\lambda_i^f\,d\lambda_i^g\).
For smooth speeds supported in facet interiors, the coefficient of the squared deformation parameter in \(J_{-r}-J_{-r}(1_T,1_T)\) is
\[
 q_r(f,g)=\frac12\big(\|f\|_{a_r}^2+\|g\|_{a_r}^2\big)+K_r(f,g),
 \qquad \|f\|_{a_r}^2=\sum_e\int a_r f_e^2\,d\nu_e. \tag{N2}
\]
To see the factors, the two linear terms in the exact bilinear expansion each integrate \(U_{r,ij}(z+tn_{ij})=a_r(z)t+o(t)\) on the corresponding strip. Smoothing of the cross kernel gives \(K_r(f,g)\) from the product of the two boundary fluxes. The same coefficient applies if the speeds satisfy \(Cf=Cg\) and the flows acquire second-order mass matching corrections, since \(U_{r,ij}\) vanishes on the supporting planes. For such corrections, add smooth strip bumps on the edges of a spanning tree to one partition. Their Gaussian mass derivative is a nonzero diagonal rescaling of the tree incidence map, which is invertible onto the centered mass vectors. The implicit function theorem gives corrections of order \(t^2\) when \(Cf=Cg\).

Write \(\mathsf A_r\) for the self-adjoint normalized kernel \(a_r^{-1/2}K_r a_r^{-1/2}\) on the ordinary \(L^2(\nu)\) facet sum, whenever bounded. For translations \(f_{ij}^v=-v\cdot n_{ij}\) and rotations \(f_{ij}^\Omega=-(\Omega x)\cdot n_{ij}\) (\(\Omega^T=-\Omega\)), one has
\[
 K_r f^v=(a_r/r) f^v,\qquad
 K_r f^\Omega=-a_r f^\Omega. \tag{N3}
\]
For translations, the smoothed score of \(T_i+tv\) at \(x\) equals that of \(T_i\) at \(x+(t/r)v\). For rotations through \(Q_t=\exp(t\Omega)\), the score of \(Q_tT_i\) at \(x\) equals that of \(T_i\) at \(Q_t^{-1}x\). Differentiate the score difference \(T_{-r}(1_{T_i}-1_{T_j})=-U_{r,ij}\) on the facet and use its gradient \(-a_r n_{ij}\). Differentiation on the partition side integrates the boundary velocity against the conditional Gaussian kernel. These operations follow by differentiating Gaussian integrals on the moving cones, with integrable Gaussian bounds on the flat faces. In \(E\), the normalized translation and rotation eigenspaces in (N3) have dimension three each. The latter dimension follows because if all rotational speeds vanish then each root normal lies in the kernel of \(\Omega\), and the roots span \(E\).

In the active space fix \(r=1/3\). Let \(m=\nu_{12}(F_{12})\), \(\ell=4m\), and \(\theta=\arctan\sqrt2\). On \(F_{12}\) put
\[
 \xi=(1,1,-1,-1)/2,\qquad \zeta=(0,0,1,-1)/\sqrt2,\qquad x=u\xi+t\zeta.
\]
The facet has \(|t|\le\sqrt2\,u\), and \(m=2\theta/(2\pi)^{3/2}\). By the translation identity in (N3), or by applying it to an arbitrary \(v\), the scalar row for \(a_r/r\) is
\[
 K_{\rm std}(x,y)=2M_{-r}(x,y)+M_{-r}(x,(23)y)+M_{-r}(x,(24)y),
 \qquad a_r(x)=r\int K_{\rm std}(x,y)d\nu_{12}(y). \tag{N4}
\]
For example the \(1j\) and \(2j\) terms in the row have identical geometric kernels on transporting their facets back to \(F_{12}\). Their signed differences on a translation field combine to \(v_1-v_2\). The self term supplies coefficient \(2\). Take \(v_1\ne v_2\) in (N3). This gives (N4).

In (N4) rescale \(y=\sigma Y\) to absorb its Gaussian square into the surface density. If \(Y\) has the conditional standard Gaussian law on the facet, its mean is \(\sqrt{\pi/2}(\sin\theta/\theta)\xi\), by polar integration. The three transports \(I,(23),(24)\) have mixture weights \(1/2,1/4,1/4\). The last two transported mean directions project to opposite multiples of \(\zeta\) on the supporting plane. Jensen's inequality therefore gives
\[
 a_r(x)\ge a_0^{\rm sc}\exp(-\beta|x|^2-\gamma u),\qquad
 a_0^{\rm sc}=r\ell/\sigma,\quad \beta=\frac1{16},\quad
 \gamma=\frac{\sqrt\pi\sin\theta}{8\theta}<\frac15. \tag{N5}
\]

## The active-space gap

For \(x\in F_e,y\in F_h\), write \(I=d_e\cdot d_h\) and
\[
 W(x,y)=\frac12\big(I\,x\cdot y-(d_e\cdot y)(x\cdot d_h)\big).
\]
This is the inner product of the skew tensor parts of \(d_e\otimes x\) and \(d_h\otimes y\). Let
\[
 p(u,t)=1-\frac{19u}{100}+\frac{3u^2}{40}+\frac{7t^2}{200}
                -\frac{u^3}{100}-\frac{ut^2}{200},\qquad
 s_0=\frac{87}{80},\quad s_1=-\frac{49}{120}.
\]
Each facet uses the transported \(u,t\), with either orientation of \(t\).
Subtract from the kernel of \(\mathsf A_r\) the operator \(\mathsf P\) with kernel
\[
 \frac{\sigma^{-3} e^{-\beta(|x|^2+|y|^2)}}{\sqrt{a_r(x)a_r(y)}}\,
          (s_0 I+s_1W(x,y))p(x)p(y).
\]
The \(I\) term factors through \(E\) and the \(W\) term through the skew tensors on \(E\), so \(\operatorname{rank}\mathsf P\le6\). Inequality (N5) ensures their factor functions are square integrable. The [interval checker](check_tetrahedral_negative_paired.py) and the integral bounds below give
\[
 \|\mathsf A_r-\mathsf P\|_{\mathrm{HS}}<33/50,\qquad
 1-\frac{m^2}{\int a_r\,d\nu_{12}\int a_r^{-1}d\nu_{12}}<143/1000. \tag{N6}
\]
In particular the normalized operator is compact. On the orthogonal complement of the six exact eigenvectors (N3) its norm is at most \(33/50\). Indeed, if there were seven independent singular vectors with singular values exceeding \(33/50\), their span would contain a unit vector in \(\ker\mathsf P\). That would contradict the bound on \(\mathsf A_r-\mathsf P\).

Use \(f_+=(f+g)/\sqrt2\), \(f_-=(f-g)/\sqrt2\), and set \(z_\pm=\sqrt {a_r}\,f_\pm\). Then \(2q_r(f,g)\) is the sum of the forms of \(1+\mathsf A_r\) at \(z_+\) and \(1-\mathsf A_r\) at \(z_-\). The first has kernel the rotations and gap at least \(17/50\) on the orthogonal complement. The constraint is \(Cf_-=0\). In ordinary facet \(L^2\), the vectors defining \(Cf_-=0\) are \(a_r^{-1/2}(d_e\cdot v)\), \(v\in E\). The translation eigenspace is spanned by \(a_r^{1/2}(d_e\cdot v)\). The three Gram matrices for these two spans and their pairing are, respectively, \(4(\int a_r^{-1}d\nu_{12})I_E\), \(4(\int a_r d\nu_{12})I_E\), and \(4m I_E\). Thus the largest squared norm of the orthogonal projection of a balanced unit vector onto the translation eigenspace equals the second expression in (N6). On that constrained subspace,
\[
 \langle z_-,\mathsf A_r z_-\rangle
 \le \left[\frac{33}{50}+\left(3-\frac{33}{50}\right)\frac{143}{1000}\right]\|z_-\|^2
 =\frac{49731}{50000}\|z_-\|^2. \tag{N7}
\]
The rotation eigenvalue \(-1\) obeys the upper bound on the nontranslation complement.

In extra coordinates of dimension \(d-3\), expand both speed fields in an orthonormal Gaussian Hermite basis, using the product facet measure and kernel. A transverse degree \(j\) replaces \(\mathsf A_r\) by \((-r)^j\mathsf A_r\) and leaves the norm unchanged. Only \(j=0\) has the mass constraint. On the active space \(-1\le\mathsf A_r\le3\). For the symmetric form at \(j=1\), equality in \(1-r\mathsf A_r\ge0\) uses the translation eigenspace, exactly the speeds of rotations mixing the active and unused coordinates. Off this eigenspace the gap is at least \(1-r(33/50)\). The antisymmetric form at \(j=1\) has gap at least \(1-r\), and both forms for degrees \(j\ge2\) have gap at least \(1-3r^2\). Consequently (N2) at \(r=1/3\) has precisely the common-rotation kernel and obeys
\[
 q_{1/3}(f,g)\ \ge\ \frac{269}{100000}
       \operatorname{dist}_{a_{1/3}}\big((f,g),\{(h,h):h\text{ a rotational speed}\}\big)^2 \tag{N8}
\]
on \(Cf=Cg\) in every \(d\ge3\).

### Interval calculations for (N6)

The squared Hilbert–Schmidt residual is bounded using (N5). Let \(x=(u,t)\), \(y=(v,z)\) in the polar facet coordinates, with radii \(R,S\) and angles in \([-\theta,\theta]\). The upper bound is the sum over ordered edge pairs of the integrals
\[
 \frac{\sigma^{-4}}{64r^2\theta^2}
 \int RS\,e^{-h(R^2+S^2)+\gamma(u+v)}
 \big[Ie^{-\tau x\cdot y}-(s_0I+s_1W)p(x)p(y)\big]^2
                 \,dR\,dS\,d\vartheta\,d\varphi,\quad h=9/16.
 \tag{N9}
\]
The pair types are self, adjacent, disjoint, with multiplicities \(6,24,6\). Orient a representative of each as \(12,12\), \(12,13\), \(12,34\). Changing an orientation multiplies the bracket by \(-1\) when appropriate and leaves the square unchanged. In the self type \(I=2\) and \(W=x\cdot y=uv+tz\). In the adjacent type \(I=1\),
\[
 d=x\cdot y=(uz+tv)/\sqrt2+tz/2,\qquad
 W=(d-(u-t/\sqrt2)(v-z/\sqrt2))/2 .
\]
In the disjoint type \(I=0\) and \(W=-tz\). These formulas follow by writing out the four centered coordinates of \(x,y\). Sign changes in \(t,z\) preserve the integration ranges.

For the mean of \(a_r/a_0^{\rm sc}\) with respect to the normalized facet law, (N4) gives the same polar integral with exponent \(-h(R^2+S^2)\), the bracket replaced by \(2e^{-\tau d_{\rm self}}+2e^{-\tau d_{\rm adj}}\), and prefactor \((16\sigma^2\theta^2)^{-1}\). Here both adjacent terms have equal integrals by reflection of the unused labels. For the reciprocal mean, (N5) and \(u\le R\) give the one-dimensional upper bound
\(\int_0^\infty R e^{-(7/16)R^2+\gamma R}\,dR\). Integration by parts and completing the square yield the error-function expression in the checker.

All four integrations in (N9) on \(0\le R,S\le12\) use the 32-point Gauss–Legendre rule with nodes and weights enclosed in Arb intervals. If a function on an interval of length \(b\) extends holomorphically to the Bernstein ellipse of parameter \(2\) and has modulus at most \(M\) there, the Gauss error is bounded by \(8bM\,2^{-2N}\). The Chebyshev coefficient bound \(2M\,2^{-j}\) follows by extracting the Laurent coefficients on \(|w|=2\) under the substitution \((w+w^{-1})/2\). Truncate before degree \(2N\), use exactness up to \(2N-1\), and bound both integrals of the remainder by their supremum times \(b\), since the quadrature weights are positive. Telescoping the four rules gives an error at most \(32(4\cdot12^2)M\,2^{-64}\), since the angular side lengths are less than \(2\).

For the elliptic bound on each single-variable slice of (N9), use \(M=4\cdot10^{11}\). In a radial slice, \(|R|\le27/2\), \(|\operatorname{Im}R|\le9/2\). For \(a=|\operatorname{Re}R|\), the Gaussian exponent contributes at most \(h(9/2)^2-h(a^2+S^2)\), and the real dot product at most \(aS\) in modulus. For an angular slice, \(|\operatorname{Im}\vartheta|\le3/4\), so the real projection of the complex radius vector has norm \(\le(13/10)R\) and the complex Hermitian norm is \(\le(8/5)R\). The real \(u\) exponent thus uses \((13/10)\gamma R\), and \(|\operatorname{Re}(x\cdot y)|\le(13/10)RS\). With \(U,V\) bounding complex Hermitian norms, the substitution
\[
 P(U)=1+.19U+.11U^2+.015U^3,\quad
 |p(x)|\le P(U),\quad |I|\le2,\quad |W|\le2UV
\]
bounds the separated polynomial piece of the squared bracket, including \(RS\), by \(26(1+U^3)(1+V^3)P(U)^2P(V)^2\) using \(|RS|\le UV\), \(|s_0 I|\le9/4\), \(|s_1W|\le UV\), and \(|a+b|^2\le2|a|^2+2|b|^2\). In the complex angular slice take \(U=(8/5)R,V=S\). In the radial slice take \(U=a+9/2,V=S\). Absorb the real linear terms against \(t^2/16\) by completing squares and retain the factors \(e^{-t^2/2}\). For the exponential piece use \(2\tau |xy|\le\tau(|x|^2+|y|^2)\) with the corresponding real-part norm bound. The checker evaluates the resulting scalar estimates and polynomial envelopes by upper-bounding each monomial \(t^j e^{-t^2/(2b)}\) by \((3jb/8)^{j/2}\) for \(j>0\), using \(e>8/3\). For the mean integral the unsquared exponential bracket permits the same argument with exponent \(\tau\) and bound \(M=10^8\).

Outside the radial square, for the exponential part of (N9) retain \(e^{-(R^2+S^2)/6}\) after absorbing \(\gamma(R+S)\), contributing at most \(1728e^{-12^2/6}\) per edge pair, including the angular area. For the polynomial part retain \(e^{-(R^2+S^2)/4}\). Absorb \(\gamma(R+S)\) into \( (R^2+S^2)/16\), and bound the separated polynomial factors against the remaining \(e^{-(R^2+S^2)/4}\) by the same monomial envelope with \(b=2\). The integral outside the square of the retained Gaussian is at most \((8/12)e^{-12^2/4}\) before the angular factor, by \(\int_R^\infty e^{-t^2/4}\,dt\le(2/R)e^{-R^2/4}\). The mean integral has exterior bound \((512/9)e^{-3\cdot12^2/8}\), obtained by keeping \(e^{-3(R^2+S^2)/8}\) and its radial factors. These are absolute errors before multiplication by the integral prefactors. The [internal replay receipt](tetrahedral_negative_paired_verification.json) records the interval sums and errors. The checker evaluates the integral bounds derived above.

## Measurable local comparison

The two active gaps in (N6)–(N7) verify the hypotheses of the [weighted gap criterion](NEGATIVE_SIMPLEX_INTERFACE_LIFT.md#nearby-measurable-pairs-in-every-ambient-dimension) at \(s=1/3<1/2\). Its normal flux bound and transverse hypercontractivity estimate apply to measurable ordinary and fractional assignments. The criterion gives (N1) with constants independent of \(d\), and equality in the minimum value exactly for common tetrahedral assignments. It uses the corresponding-mass constraints stated above. The positive same-partition functional has a different second variation.

Global minimization at four cells and the [large-alphabet negative voting interval](NEGATIVE_VOTING_INTERVAL.md) remain separate questions. The proof and interval replay have internal review. Independent mathematical review and formal verification remain open.
