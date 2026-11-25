# The period obstruction for stable leading sphere phases

Let \(K\) have characteristic zero, and let
\[
 Q=\{t^2+2xy=1\},\qquad
 \Lambda(f)=\operatorname{CT}_x\frac12\int_{-1}^{1}
 f\bigl(x,(1-t^2)/(2x),t\bigr)\,dt.
\]
Sphere degree is the least degree of a polynomial representative on \(Q\).
For a phase \(f\) of degree \(d>0\), let \(H_d\) be its leading symbol on
the homogeneous cone and \(h_d(u,v)=H_d(u^2,-v^2/2,uv)\) its binary form.

**Theorem.** If every geometric root of \(h_d\) has multiplicity strictly
less than \(d\), then \(\Lambda(f^m)\ne0\) for infinitely many \(m\).

Together with the [exact-half-root obstruction](SPHERICAL_LEADING_BINARY_FORM_MOMENT_OBSTRUCTIONS.md#4-proof-of-the-elementary-obstructions),
this shows that every nonconstant eventually pure sphere phase has a leading
binary form with a root of multiplicity greater than its sphere degree.

## Compactification and extension of the volume form

The coefficient field of \(f\) can be taken finitely generated over
\(\mathbb Q\) and embedded in \(\mathbb C\). Rational moment identities and
geometric root multiplicities persist under this embedding. Let
\(\bar Q\subset\mathbb P^3\) be the projective quadric with boundary conic
\(D\) at infinity. The affine volume form
\[
 \omega=\operatorname{Res}_{Q}\frac{dx\wedge dy\wedge dt}{t^2+2xy-1}
\]
is nowhere vanishing on \(Q\). Adjunction gives \(K_{\bar Q}\sim-2D\), and
the divisor of \(\omega\) is supported on the irreducible conic \(D\), so
\(\operatorname{div}(\omega)=-2D\). In particular near a boundary point,
local coordinates \((u,z)\) with \(D=\{u=0\}\) give
\(\omega=u^{-2}g(u,z)\,du\wedge dz\), where \(g\) is a holomorphic unit.

The rational pencil defined by \(f\) has local expression
\(F(u,z)/u^d\). On \(u=0\), the local function \(F\) represents the leading
symbol as a section of \(O_{\bar Q}(d)|_D\). At a root of multiplicity
\(r>0\), Weierstrass preparation writes \(F\) as a holomorphic unit times a
monic degree-\(r\) polynomial in \(z\) on a smaller bidisc. For fixed \(u\),
every point of \(\{|F(u,z)|\le\delta\}\) lies within distance
\(C\delta^{1/r}\) of some root of that polynomial. The planar area of the
sublevel set is therefore at most \(C'\delta^{2/r}\), uniformly in \(u\).
For any bounded \(M>0\), this gives finite volume for \(\omega\) near that
boundary point on the region \(|f|\le M\), since
\[
 \int_{|u|<\epsilon}|u|^{-4}\,
       \operatorname{Area}\{z:|F(u,z)|\le M|u|^d\}\,d^2u
 \le C_M\int_{|u|<\epsilon}|u|^{-4+2d/r}\,d^2u<\infty .
\]
Here finite volume means integrability of the positive density determined
by \(\omega\wedge\bar\omega\). The last inequality uses \(r<d\). At points
where the leading section does not vanish, the bounded-\(f\) region stays
away from \(D\) in a sufficiently small neighborhood.

Resolve the base ideal of the pencil on the smooth projective surface, by
blowups supported on the boundary. This yields a smooth projective surface
\(\widetilde Q\) on which \(f\) extends to a morphism to \(\mathbb P^1\).
Set \(Y=\widetilde Q\setminus f^{-1}(\infty)\), so
\(f:Y\to\mathbb A^1\) is proper and nonconstant and the affine quadric still
embeds in \(Y\). In a relatively compact chart of \(Y\) the function \(f\)
is bounded. The volume estimate and the isomorphism away from the boundary
give finite squared integral for the pulled-back form off that boundary.
A genuine pole of a meromorphic top form along a divisor has divergent
squared integral near its generic smooth points, as follows by integrating
\(|u|^{-2k}\) in a transverse complex disk for a pole of order \(k\ge1\).
Thus the pulled-back rational form has no divisorial poles on \(Y\). Since
\(Y\) is smooth, it is a global algebraic holomorphic top form there.

## Applying the proper period obstruction

The [proper holomorphic seed lemma](TWO_PAIR_IMAGE_CUBIC_STRATIFICATION.md#a-proper-holomorphic-volume-seed-cannot-have-constant-exponential-period)
applies to this \(f\) and \(\omega\). The same proof rules out any nonzero
*polynomial* exponential period. To see this, in its notation let \(m\) be
the direct-image class of \(\omega\). It belongs to \(V^{>0}M\) by negativity
of the Bernstein--Sato roots and proper direct image of the graph module.
The cyclic Fourier module of \(m\) acts on its period by the Weyl operators.
If that period is a nonzero polynomial \(P(\tau)\), differentiating \(P\) to
its top coefficient shows that its Weyl submodule is all of
\(\mathcal O_{\mathbb A^1_\tau}\). Inverse Fourier transform would then give
a quotient of the cyclic submodule of \(m\) onto \(\delta_0\), with nonzero
image of \(m\). The induced filtration on submodules and functoriality of
the \(V\)-filtration prohibit this, since \(V^{>0}\delta_0=0\).

The period evaluation in this argument uses the Fourier-transformed
relative complex. The graph action \(\partial_t\mapsto\tau\) changes the
relative de Rham differential to \(d_Y-\tau\,df\wedge\) and
\(\partial_\tau\) acts by \(\partial_\tau-f\). Specialization at each \(\tau\)
maps hypercohomology relations to relations for this specialized complex.
Multiplication by \(e^{-\tau f}\) on the analytic complex maps those to
ordinary de Rham relations. Thus integration on a compact closed cycle
respects the Weyl relations, including when the algebraic relation is
represented using Čech cochains rather than a global primitive. The
filtration inputs are Chen--Dirks--Mustaţă,
[*An introduction to \(V\)-filtrations*](https://arxiv.org/abs/2402.07816v1),
Example 5.19, Theorem 8.1, Proposition 3.19 and Corollary 3.20, in the
decreasing convention of the proper seed lemma.

On the real sphere \(\Gamma=\{y=\bar x,\ 2|x|^2+t^2=1,\ t\in\mathbb R\}\),
the integral of \(\omega\) is a nonzero constant and normalized integration
is exactly \(\Lambda\), by the
[sphere residue formula](GLOBAL_SPHERE_PURITY_RECOGNITION.md#1-a-closed-rational-period).
The entire function \(\int_\Gamma e^{-\tau f}\omega\) has coefficients
\((\int_\Gamma\omega)(-1)^m\Lambda(f^m)/m!\) by uniform convergence of the
exponential expansion on this compact cycle. If all sufficiently large
moments vanished, it would be a polynomial with nonzero constant term,
contradicting the proper period obstruction. This proves the theorem.
