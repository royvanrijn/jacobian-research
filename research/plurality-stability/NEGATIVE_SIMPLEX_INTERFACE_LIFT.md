# Weighted interface form for opposite simplex partitions

Fix \(k\ge3\), write \(n=k-1\), and realize the regular simplex fan
\(T=(T_i)_{i=1}^k\) in \(E=\mathbf1^\perp\subset\mathbb R^k\). Let
\(F_{ij}\) be the common facet of \(T_i,T_j\),
\(n_{ij}=(e_i-e_j)/\sqrt2\) point into \(T_i\), and let
\(d\nu_{ij}=\phi_n\,d\mathcal H^{n-1}|_{F_{ij}}\). For \(0<s<1\)
consider the bilinear Gaussian stability

\[
 P_s(A,B)=\sum_i\langle 1_{A_i},\,T_{-s}1_{B_i}\rangle_{L^2(\gamma)}
\]

over pairs having \(\gamma(A_i)=\gamma(B_i)\) for each \(i\), with
the common mass vector allowed to vary. Reflecting \(B\) through the
origin converts this to the positive-correlation bilinear problem with
reference pair \((T,-T)\). This is the correspondence in
[Heilman, Problem 7.1 and Lemma 7.6](https://arxiv.org/html/2209.11216#S7).

## Weighted Hessian and compactness

Orient speeds by the change in the indicator of the first label,
so \(\delta1_{T_i}=\lambda_i^f=\sum_{j\ne i}f_{ij}\nu_{ij}\) with
\(f_{ji}=-f_{ij}\). Put

\[
\begin{aligned}
 a_s(x)&=-\partial_{n_{ij}}T_{-s}(1_{T_i}-1_{T_j})(x),
    &&x\in F_{ij},\\
 D_s(f,g)&=\sum_{i<j}\int_{F_{ij}}a_s(x)f_{ij}(x)g_{ij}(x)\,d\nu_{ij}(x),\\
 K_s(f,g)&=\sum_i\iint M_{-s}(x,y)\,d\lambda_i^f(x)\,d\lambda_i^g(y),\\
 (Cf)_i&=\sum_{j\ne i}\int_{F_{ij}}f_{ij}\,d\nu_{ij}.
\end{aligned}
\tag{NL1}
\]

Here \(M_{-s}\) is the Mehler density relative to the product Gaussian
measure. Reflection in the pair plane interchanges \(T_i,T_j\) and makes
the smoothed score difference zero on the whole plane. For a displacement
out of \(T_i\), its derivative is \(a_s\). Integration over a thin normal
strip therefore contributes \(D_s(f,f)/2\) to the coefficient of the
squared displacement parameter. The mixed term is \(K_s(f,g)\).
The second derivative at the common simplex \(T\) is therefore

\[
 D_s(f,f)+D_s(g,g)+2K_s(f,g)
   =(D_s+K_s)(u,u)+(D_s-K_s)(v,v),\qquad
 u=(f+g)/\sqrt2,\quad v=(f-g)/\sqrt2 .
 \tag{NL2}
\]

The matched-mass constraint is \(Cv=0\).
The strip argument applies first to compactly supported smooth speeds
inside facets. Small additional strips on edges of a spanning tree of
the labels correct the higher-order mass error: their mass derivatives
are nonzero multiples of the tree incidence vectors, which span
\(\mathbf1^\perp\), so the implicit function theorem applies. The
original linear variation of the objective vanishes because the
reference scores tie on every pair plane.

Let \(\sigma=\sqrt{1-s^2}\), \(q=z\cdot n_{ij}\), and
\(\beta=s^2/(2\sigma^2)\). Differentiation of the reflected score gives

\[
 a_s(x)=\frac{2s}{\sigma}\int_{T_i}q\,\phi_n(z+sx/\sigma)\,dz.
 \tag{NL3}
\]

Enlarging to \(q\ge0\) gives \(a_s(x)\le 2s\phi_1(0)/\sigma\).
Integrating over a fixed compact ball in the interior of \(T_i\) scaled
by \((1+|x|)^{-1}\) gives, for some \(c>0\),

\[
 c(1+|x|)^{-(n+1)}e^{-\beta|x|^2}\le a_s(x).
 \tag{NL4}
\]

Indeed \(q\) is at least a positive constant times \((1+|x|)^{-1}\)
on that ball, its volume scales as \((1+|x|)^{-n}\), and the cross term
in \(|z+sx/\sigma|^2\) is bounded there. Along every interior facet ray
\(x=t x_0\) with \(x_0\ne0\), \(a_s(x)\to0\). The shift \(-sx/\sigma\)
then goes to infinity at linearly growing distance from \(T_i\), and its
normal component is zero. In (NL3) the event is contained in the Gaussian
tail outside that growing distance, with \(q\) equal to the centered
normal coordinate.

From now on take \(0<s<1/2\). The Hilbert space of speeds has norm
\(\|f\|_{D_s}^2=D_s(f,f)\). The bound (NL4) shows
\(\int a_s^{-1}d\nu_{ij}<\infty\), so \(C\) is continuous in this norm.
The interaction is compact on the weighted space. In fact after the
unitary multiplication \(f\mapsto\sqrt{a_s}f\) its kernel on any pair
of facets is the incidence coefficient times
\(M_{-s}(x,y)/\sqrt{a_s(x)a_s(y)}\). The square of this kernel times
the two Gaussian surface densities is bounded by a polynomial factor
times

\[
 \exp\!\left[-\frac{1-2s}{2(1-s^2)}(|x|^2+|y|^2)\right].
 \tag{NL5}
\]

This follows from \(2|x\cdot y|\le |x|^2+|y|^2\) and (NL4).
Each facet block is Hilbert--Schmidt, and their finite sum is
self-adjoint. Compactness and (NL1) extend the forms by density to the
whole weighted space. The same compactness implies a positive weighted
gap whenever a nonnegative block \(D_s\pm K_s\) has a specified kernel,
also on the closed subspace \(\ker C\). In normalized coordinates it
is the identity plus or minus a compact self-adjoint operator.

## Dimension lift for the paired form

Let \(\mathcal T\) be the active translation space
\(t_{ij}^w=w\cdot n_{ij}\) for \(w\in E\), and let \(\mathcal R\) be the
active rotation space \(r_{ij}^{\Omega}(x)=-(\Omega x)\cdot n_{ij}\)
for skew-adjoint \(\Omega:E\to E\). Translation and rotation satisfy the
mixed form identities

\[
 K_s(t,h)=s^{-1}D_s(t,h),\qquad
 K_s(r,h)=-D_s(r,h)
 \tag{NL6}
\]

for every \(h\). For translation use
\(\nabla1_{T_i}=\sum_{j\ne i}n_{ij}\mathcal H^{n-1}|_{F_{ij}}\) and
\(\partial_w T_{-s}1_{T_i}=-s T_{-s}(\partial_w1_{T_i})\), with the
Gaussian surface convolution in (NL1). The gradient of the smoothed
pair difference on its pair plane is \(-a_s n_{ij}\).
For rotation, differentiating the indicators of the rotated fan gives
\(\lambda_i^r\), and commutation with rotation gives the smoothed
derivative \(-(\Omega x)\cdot\nabla T_{-s}(1_{T_i}-1_{T_j})\)
on \(F_{ij}\), which equals \(-a_s r_{ij}^\Omega(x)\).
Gaussian integration on the facets justifies these derivatives,
and (NL5) extends the identities from compact tests.

**Dimension lift.** Suppose the active forms \(D_s+K_s\) on all speeds
and \(D_s-K_s\) on \(\ker C\) have weighted gaps \(c_+>0\) modulo
\(\mathcal R\) and \(c_->0\) with zero kernel, respectively. Then the
paired form (NL2) on the cylindrical simplex in every \(d\ge k-1\)
has gap at least \(\min\{c_+,c_-,1-s\}\) in its weighted product norm,
modulo the common rotations. Its kernel is precisely those common
rotations.

To prove this, split any active \(h=t+b\) with \(t\in\mathcal T\) and
\(Cb=0\). The integral of every facet is the same positive number
\(\nu_F\), and \(Ct^w=k\nu_Fw/\sqrt2\), so this split exists and is
unique. By (NL6) all mixed terms with \(t\) vanish in \(D_s-sK_s\)
and

\[
 (D_s-sK_s)(h,h)
   =(1-s)D_s(b,b)+s(D_s-K_s)(b,b)
   \ge(1-s)\operatorname{dist}_{D_s}(h,\mathcal T)^2 .
 \tag{NL7}
\]

Together with active \(D_s+K_s\ge0\), this puts the spectrum of the
weighted interaction operator in \([-1,s^{-1}]\).
For unused coordinates \(z\in\mathbb R^{d-n}\), the local coefficient
is independent of \(z\). Expanding in Gaussian Hermites of \(z\)
multiplies the interaction in degree \(m\) by \((-s)^m\).
Only the constant degree carries the mass condition.
In degree one the symmetric block is \(D_s-sK_s\), with the gap
and translation kernel in (NL7). The antisymmetric block is
\(D_s+sK_s\ge(1-s)D_s\). In every degree \(m\ge2\),
\(|(-s)^m K_s(h,h)|\le s^{m-1}D_s(h,h)\) gives the common
lower bound \(1-s\) in both blocks. The degree-one translation speeds
times linear \(z\) are exactly the rotations mixing the active and
unused spaces. Together with degree-zero active rotations they
exhaust the common rotational kernel. The tensor Hermite decomposition
and boundedness of the operators justify summing all blocks.

## The complete rotation hook at negative correlation

**Active rotation hook.** For every \(k\ge3\) and \(0<s<1/2\),
\(D_s+K_s\) on the \(S_k\) isotypic space of \(\Lambda^2 E\) has a
strictly positive weighted gap modulo \(\mathcal R\), and \(D_s-K_s\)
has a strictly positive weighted gap there.

Use the full rotation-hook multiplicity model and coordinates from
[(S52)--(S53)](SIMPLEX_INTERFACE_SYMMETRY.md#the-complete-hook-kernel).
They use only the Gaussian kernel and the geometry and apply with
\(\tau=-s/(1-s^2)<0\). In detail, on the chamber
\(W=\{x\in F_{12}:x_2>x_3>\cdots>x_k\}\) use gap coordinates
\(\eta_j=(v_j-v_{j+1})/k\), \(2\le j\le k-1\). Their inverse \(V\)
has columns \(r_m\) with entries \(k-m\) in the first \(m\) positions
and \(-m\) in the others. The metric \(J=V^TV\) has positive entries.
For an input column \(r_m\) and permutation \(\pi\) of the untied labels,
put \(z=\pi y\), \(v=\pi r_m\), \(b=v_2=k-m\),
\(E_0=e^{\tau x\cdot z}\), and \(E_l=e^{\tau x\cdot(2l)z}\).
The output is the sum over \(\pi\) of the gap coordinates of
\(2E_0v+\sum_{l\ge3}E_l(I-2w_lw_l^T)v\), where
\(w_l=(e_1+e_2)/2-e_l\). A common positive Gaussian factor
multiplies the sum.

Every entry of this folded hook kernel is strictly negative.
For its \(j=2\) row the self and \(l=3\) terms combine as
\(2(b-v_3)(E_0-E_3)\le0\), since
\(x\cdot z-x\cdot(23)z=(x_2-x_3)(z_2-z_3)>0\).
For \(l\ge4\), pairing with \((3l)\pi\) gives coefficient
\(v_l-v_3\) times an exponential difference with exponent difference
\(\tau(x_2-x_3)(z_l-z_3)\). This contribution is nonpositive since
\(v_l-v_3\) has the sign of \(z_l-z_3\) when nonzero.
Some \(\pi\) has \(v_3=-m\), making the first inequality strict.
For \(j\ge3\), pair with \(p\pi\), \(p=(j,j+1)\). The self and
nonexceptional terms have coefficient \(v_j-v_{j+1}\), with a factor
of two for the self term. Their differences have exponent difference
\(\tau(x_j-x_{j+1})(z_j-z_{j+1})\) and give nonpositive contributions.
The exceptional terms combine as

\[
 (2b-v_j-v_{j+1})[E_j(z)-E_{j+1}(pz)+E_j(pz)-E_{j+1}(z)].
\]

The two exponent differences in brackets are
\(\tau(x_j-x_{j+1})(z_2-z_{j+1})\) and
\(\tau(x_j-x_{j+1})(z_2-z_j)\), both negative.
For \(m\ge3\) the self pairing is strict for a permutation giving
\(v_j\ne v_{j+1}\), and for \(m=2\) the exceptional coefficient is
strictly positive. This proves the strict sign in all rows.

Let \(A\) be minus this folded kernel normalized by \(a_s^{-1/2}\)
at each argument. It is compact and self-adjoint in the \(J\)
inner product by (NL5). The rotation field in multiplicity coordinates
is \(v(x)=x\). By (NL6), \(\psi=\sqrt{a_s}\,V^{-1}x\) is a right
eigenfield of \(A\) with eigenvalue \(1\), and \(\chi=J\psi\) is a
left eigenfield in ordinary coordinates. They have strictly positive
components on \(W\) and belong to \(L^2\). For a real eigenfield \(u\)
with eigenvalue \(\lambda\), the triangle inequality gives

\[
 |\lambda|\int_W\chi\cdot|u|\,d\nu_{12}
 \le\int_W\chi\cdot A|u|\,d\nu_{12}
 =\int_W\chi\cdot|u|\,d\nu_{12}.
\]

Equality forces all nonzero components of \(u\) to have the same
sign almost everywhere, by strict positivity of every kernel entry.
It excludes eigenvalue \(-1\) and makes the eigenspace at \(1\)
one-dimensional: a second independent eigenfield could be chosen
orthogonal to \(\psi\) in the \(J\) metric. Compactness then gives a
spectral gap below \(1\) on its orthogonal complement and above \(-1\)
on the entire space. Restoring the representation coordinate makes
the one ground field exactly the geometric active rotation copy,
giving both weighted gaps.

The dimension lift requires the whole active pair of forms to satisfy
its hypotheses. The hook proof supplies the rotation block. At
positive correlation with one partition the interaction has its
positive Mehler coefficients added at both parities and is subtracted
from a different local coefficient. With two partitions at negative
correlation the two forms in (NL2) use opposite signs of the
alternating Hermite series.

## Nearby measurable pairs in every ambient dimension

**Weighted gap criterion.** Assume the two active gaps in the dimension
lift at \(0<s<1/2\). There exist \(c,\eta,\varepsilon>0\), depending
only on \(k,s\) and the gaps, such that every \(d\ge n\),
\(|r-s|<\varepsilon\), and pair of measurable simplex-valued maps
\(A,B\) with \(\int A_i\,d\gamma_d=\int B_i\,d\gamma_d\) satisfy

\[
 P_r(A,B)-P_r(T,T)\ge c\,\operatorname{dist}_1((A,B),\mathcal O)^2
 \quad\text{if }\operatorname{dist}_1((A,B),\mathcal O)<\eta .
 \tag{NL8}
\]

The functional uses \(A_i,B_i\) in place of the indicators above.
Here \(\mathcal O\) is the orbit of the common cylindrical simplex
under simultaneous orthogonal transformations and label permutations.
Distance is the infimum of summed Gaussian \(L^1\) differences over
both maps and all labels. In this neighborhood equality in the
minimum value occurs exactly on \(\mathcal O\).

**Moment alignment.** Choose the reference simplex \(T\) to minimize
the Frobenius distance of its moment matrix from the average of the
moment matrices of \(A,B\). Each row of a moment matrix is
\(\int (x,z) A_i(x,z)\,d\gamma_d\). This minimum exists on the
compact orbit in each dimension. If the original pair approaches
\(\mathcal O\) in \(L^1\), then it also approaches this chosen \(T\).
In fact the squared Frobenius error in the moment matrices is
bounded by the \(L^2\) squared error of the maps by Hermite projection,
which is at most their summed \(L^1\) error. Thus the minimizing
reference moment matrix approaches that of any initially aligning
simplex. The normalized pair differences of the rows give its facet
normals. Two unit normals whose Euclidean distance tends to zero
give halfspaces whose Gaussian symmetric difference is their angle
divided by \(\pi\). The union bound over pair planes proves the
alignment claim with bounds independent of \(d\).

Write \(\mathcal D f\) for the matrix of first moments of the facet
measures in (NL1), including the unused coordinates. Differentiating
the nearest-orbit condition shows that the average moment error is
orthogonal to the tangent space of the regular moment orbit. This
constraint provides a uniform coercive gauge for the symmetric speed.
Indeed \(C\) and \(\mathcal D\) are bounded in the weighted norm
uniformly in \(d\): for the moment bound test each row against a unit
vector in \(E\oplus E^\perp\), use Cauchy--Schwarz with the weight
\(a_s\), and integrate its squared scalar product against
\(a_s^{-1}d\nu\,d\gamma_{d-n}\). The active contribution is finite by
(NL4), and the unused contribution is a standard one-dimensional
Gaussian second moment. There are only \(k\) rows.

In active and unused coordinates the regular moment matrix has the
form \(a_*[I_E,0]\) in orthonormal coordinates on the centered label
space, with \(a_*>0\) depending only on \(k\).
On common rotations, the moment map sends the speed of a skew matrix
\(\Omega\) to \(a_*[I_E,0]\Omega^T\), by differentiation of the rotated
moments. The active and mixed blocks of \(\Omega\) determine the speed,
and their Frobenius norms bound its weighted norm using (NL3) and
Gaussian second moments on the facets. The restriction of
\(\mathcal D\) to rotations therefore has an inverse from the tangent
space with norm bounded independently of \(d\).
Consequently, if \(\mathcal D u\) has zero tangent projection and
\(u=u_0+u_{\rm rot}\) is its weighted orthogonal decomposition into
nonrotational and rotational parts, the tangent projection of
\(\mathcal D u_0\) bounds \(\|u_{\rm rot}\|_{D_s}\) by a fixed
multiple of \(\|u_0\|_{D_s}\).
The dimension lift gives some \(\delta_*>0\), independent of \(d\),
such that on this symmetric gauge and on the antisymmetric mass gauge
the paired form satisfies

\[
 D_s(f,f)+D_s(g,g)+2K_s(f,g)
 \ge\delta_*[D_s(f,f)+D_s(g,g)] .
 \tag{NL9}
\]

**Control of discarded transfers.** Write \(e^A=A-T,e^B=B-T\),
\(h^A_{ij}=1_{T_i}A_j\), and similarly for \(B\). With
\(V_{r,ij}=-T_{-r}(1_{T_i}-1_{T_j})\), set
\(L=\sum_{i\ne j}\int(h^A_{ij}+h^B_{ij})V_{r,ij}\,d\gamma_d\).
Direct expansion gives the bilinear deficit as \(L+\sum_i
\langle e^A_i,T_{-r} e^B_i\rangle\). Reflection shows \(V_{r,ij}>0\)
on the interior of \(T_i\). Choose \(b>s\) with \(b<1/2\) and choose
exponents \(a>b/[2(1+b)]\), \(\beta'>b^2/[2(1-b^2)]\) with
\(2a+\beta'<1/2\), possible because \(2b<1\). Take \(r\) in a
compact neighborhood of \(s\) within \((0,b)\). On \(T_i\),
\(t=x\cdot n_{ij}\ge0\), and the reflected formula (N4) in the
[small negative correlation argument](NEGATIVE_SMALL_CORRELATION.md#exact-negative-correlation-scores)
yields

\[
 V_{r,ij}(x)\ge c_0 e^{-\beta'|x|^2}\min(t,1).
 \tag{NL10}
\]

For this bound integrate over the scaled interior ball used in (NL4).
The reflected density is bounded below there by a constant times
\(e^{-r^2|x|^2/(2(1-r^2))}\). The factor
\(1-e^{-2rtq/\sqrt{1-r^2}}\) on that ball bounds a constant times
\(\min(t,1)/(1+|x|)\), and its volume supplies another fixed power
of \(1+|x|\) in the denominator. The slack in \(\beta'\) absorbs
these powers uniformly in \(r\).

Use \(W(x)=(1+|x|^2)e^{a|x|^2}\). For an ordered transfer \(h\),
put \(J_h(z)=\int W(x)h(x,z)d\gamma_n(x)\). The slice mass-cost
argument in [(N7)](NEGATIVE_SMALL_CORRELATION.md#transfer-estimates-in-the-active-coordinates)
gives \(\|J_h\|_2\le C\sqrt L\) and
\(\|J_{h1_{\{|x|>R\}}}\|_2\le\omega(R)\sqrt L\) with
\(\omega(R)\to0\). In that argument the integrable weight is now
\(W^2 e^{\beta'|x|^2}\), which is Gaussian integrable by the choice
of \(a,\beta'\). On \(|x|\le R\), the part with \(t>\xi>0\) has
\(\|J_h\|_1=O_{R,\xi}(L)\) by (NL10). Interpolation with its
\(O(\sqrt L)\) bound in \(L^2(d\gamma_{d-n})\) gives
\(o_{R,\xi}(\sqrt L)\) in \(L^p\), where \(p=1+b<2\), as \(L\to0\).
Finally, on the part with \(t\le\xi\) whose normal projection \(y\)
lies off \(F_{ij}\), some untied height lies between \(y_i\) and
\(y_i+\xi/\sqrt2\). In the projected compact ball this is a union of
strips whose plane measure tends to zero with \(\xi\). The same
normal mass-cost bound and tangential Cauchy--Schwarz with this
restricted measure give \(o_{R,\xi\downarrow0}(\sqrt L)\) in \(L^2\).

These estimates control the interactions in every unused dimension.
The active Mehler density obeys \(|M_{-r}(x,y)|\le C W(x)W(y)\).
In fact its exponent is at most
\(r(|x|^2+|y|^2)/[2(1+r)]\). After integration in the active
coordinates, the absolute interaction of two transfers is bounded
by \(C\langle J_h,T_{-r}J_{\tilde h}\rangle\le C\|J_h\|_p
\|J_{\tilde h}\|_p\). The last step uses
\(T_{-r}=T_{\sqrt r}T_{\sqrt r}\mathcal I\), where
\(\mathcal I\) reflects the unused coordinates, and Gaussian
hypercontractivity from \(L^{1+b}\) to \(L^2\), since \(r<b\).
This is Nelson's theorem in the form stated as Theorem 1.1 in
[Gross, *Hypercontractivity over complex manifolds*](https://doi.org/10.1007/BF02392573).
The same discarded-transfer bounds pass masses and first moments.
For unused moments, test against a unit vector \(w\) there and
apply Hölder with the \(L^{p/(p-1)}\) norm of \(w\cdot z\), which is
dimension independent.

**Passage to the facet speeds.** Consider any sequence with
\(r\to s\) and the pair approaching \(\mathcal O\). Dimensions can
vary with the sequence. Work in the moment alignment above, and
consider terms with \(L>0\). The scores are bounded and the
transfers tend to zero in \(L^1\), so \(L\to0\).
Keep only the transfers in \(|x|\le R,\ 0\le t\le\xi,\ y\in F_{ij}\).
Project them normally to \(y\), preserving \(z\), and divide by
\(\sqrt L\). Relative to \(d\nu_{ij}(y)d\gamma_{d-n}(z)\) their
densities are

\[
 \theta_{ij}(y,z)=\frac{1}{\phi_1(0)\sqrt L}
   \int h_{ij}^{\rm kept}(y+t n_{ij},z)\phi_1(t)\,dt .
 \tag{NL11}
\]

Take the difference of the opposite ordered densities for each
pair and call the resulting two speeds \(f,g\). On the retained compact
sets the score divided by \(t\) approaches \(a_s(y)\) uniformly as
\(r\to s,\ \xi\to0\), by differentiating the reflected scores at
\(t=0\). The weight is bounded below there. The inequality that a
density bounded by \(\phi_1(0)\) on \(t\ge0\) has first moment at
least its squared mass divided by \(2\phi_1(0)\) gives
\(D_s(f,f)+D_s(g,g)\le2+o(1)\). Here and below each error can be
made arbitrarily small by choosing \(R\) large, then \(\xi\) small,
then going far enough along the sequence. The inequality uses
\((\theta_{ij}-\theta_{ji})^2\le\theta_{ij}^2+\theta_{ji}^2\).

The retained interaction divided by \(L\) changes by
\(O_R(\xi)\) upon normal projection, by uniform differentiability of
the active kernel on compact sets and the \(L^2\) slice bounds above.
Its difference from \(K_s(f,g)\) also tends to zero for fixed \(R,\xi\).
In the active variables the Mehler densities converge uniformly on
the retained supports. In the unused variables the \(L^2\) operator
norm difference between \(T_{-r}\) and \(T_{-s}\) is bounded by
\(\sup_{m\ge0}|r^m-s^m|=O_b(|r-s|)\), by Hermite expansion. The
retained fluxes have bounded \(L^2\) norm on the compact planes
since \(a_s\) is bounded below there. Cauchy--Schwarz in their finite
surface measures and in \(z\) justifies applying that operator
difference. The discarded transfer estimates therefore give
\(\sum_i\langle e^A_i,T_{-r}e^B_i\rangle/L=K_s(f,g)+o(1)\).

The original pair has equal masses, and the original average moments
obey the tangent gauge. Projection preserves mass of each retained
transfer and its unused first moments, changing its active moment
by at most \(O(\xi)\) times the normalized mass. The discarded
transfer bounds consequently give \(C(f-g)=o(1)\) and tangent
projection \(\Pi_{\rm tan}\mathcal D(f+g)=o(1)\). Subtract a
translation from the antisymmetric speed to impose the first gauge,
using the inverse of \(C|_{\mathcal T}\) from (NL7). Subtract a
common rotation from the symmetric speed to impose the second, using
the uniform moment inverse above. Both corrections are \(o(1)\)
in weighted norm. Boundedness of \(K_s\), uniformly under the dimension
lift, shows that these corrections leave the preceding estimates
unchanged up to \(o(1)\).

Apply (NL9) to the corrected speeds, decreasing \(\delta_*\) so that
\(0<\delta_*<1\). The normalized bilinear deficit is at least
\(1-(1-\delta_*)(2+o(1))/2+o(1)=\delta_*+o(1)\).
Thus on a uniform neighborhood of the orbit and correlation \(s\),
the deficit is at least \(\delta_* L/2\). A sequence of
violations with shrinking radii would contradict this lower limit.
Finally \(\operatorname{dist}_1((A,B),\mathcal O)^2\le C L\) follows
from the slice estimates on all transfers and \(W\ge1\). If \(L=0\)
all transfers vanish because the scores are strictly positive in the
cell interiors. This proves (NL8) and its equality statement.
