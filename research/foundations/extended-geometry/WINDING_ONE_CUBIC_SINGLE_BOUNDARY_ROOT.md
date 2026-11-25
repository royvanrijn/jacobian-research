# The single-root condition on the logarithmic boundary

Use the regular winding-one cubic sphere phase
\(f=x^{-1}(C+Bt+Dt^2+Et^3)\) and sphere functional \(\Lambda\) of the
[logarithmic boundary theorem](WINDING_ONE_CUBIC_LOGARITHMIC_BOUNDARY.md).
Suppose \(E(0)\ne0\), \(\deg E\le5\), and \(E\) has a root \(r\) of
multiplicity exactly four. Let \(v=x-r\) and write
\[
 E=e_4 v^4+O(v^5),\qquad
 D=d_3 v^3+O(v^4),\qquad
 B=b_2 v^2+O(v^3),\qquad
 C=c_1 v+O(v^2).
\]
Here \(r,e_4\ne0\). These orders for \(D,B,C\) are the necessary
divisibilities supplied by the logarithmic boundary theorem.

**Proposition.** Put
\[
 g(w)=e_4 w^4+d_3 w^3+b_2 w^2+c_1w.
\]
If \(g\) has only one distinct geometric root, then
\(\Lambda(f^m)\ne0\) for infinitely many positive integers \(m\).
In the displayed coordinate the sole root is necessarily zero, so the
excluded subcase is \(v^4\mid D,\ v^3\mid B,\ v^2\mid C\).

## The finite boundary curve and its residue

Work over \(\mathbb C\) after embedding a finitely generated coefficient
field containing the roots in question. Use the compactification of the
logarithmic boundary theorem. Its ordinary blowup chart at \(r\) has
\(u=1/t,\ v=uw\). In this chart
\[
 f=u\frac{g(w)+O(u)}{r+uw},\qquad
 \omega=\frac{du\wedge dw}{u(r+uw)}
\]
up to a common nonzero constant in the volume convention. The vertical
boundary curve \(V\) is \(u=0\), with coordinate \(w\) and logarithmic
residue \(dw/r\). It meets other components of the zero fibre at
finite \(w\) only over roots of \(g\). The companion toric chart
\(v=s^3z,\ u=s^4z\) has \(V=\{z=0\}\), \(w=1/s\), and its horizontal
polar divisor \(H\) is \(s=0\). There \(f=z\) times a unit near the
intersection, with \(\partial_z f\ne0\). Thus \(f\) is a transverse
coordinate on \(H\) at \(H\cap V\).

Resolve the union of the central fibre and \(V\) near the finite
part of \(V\), leaving the already resolved chart at \(H\cap V\).
This can be done by point blowups outside the affine sphere. There the
only polar component before blowing up is the simple pole along \(V\).
A blowup centred on one simple polar component adds one Jacobian factor,
so the new exceptional curve is regular for the volume. Thus the strict
transform of \(V\) is still the sole vertical polar component from this
chart. Its residue is unchanged. The other polar charts over the central
fibre, if any, are covered by the horizontal-pole lemma of the
logarithmic boundary proof. In particular any other finite root of \(E\)
has multiplicity at most one. The strict transform of \(V\) meets the
other components of the central fibre away from \(H\) only over the
original roots of \(g\).

## Correction in the filtered graph complex

Allow meromorphic poles along the boundary polar divisors of the model
and use the graph module of \(f\) with these localized coefficients.
Use \(b\) for the target coordinate and \(\delta_f\) for the graph
generator. The relative de Rham differential in the graph presentation
is \(\nabla=d-\partial_b\,df\wedge\). It preserves the filtration
along \(b=0\), because differential operators in the source directions
preserve the target filtration. Regular coefficients times \(\delta_f\)
belong to \(V^{>0}\) by Bernstein-root negativity. The same assertion
holds with horizontal meromorphic poles on the charts of the
horizontal-pole lemma.

If \(w_0\) is the sole finite zero of \(g\), take the primitive
\(h=(w-w_0)/r\) for the residue on \(V\). It has zero at every
nonhorizontal crossing of \(V\) with the central fibre, and its only
pole on \(V\) is at \(H\cap V\). More generally, this choice of
primitive works if the limiting polynomial in the finite chart is
\(e_4(w-w_0)^4\) in another affine coordinate with the same point at
infinity.

Near a finite crossing take normal-crossing coordinates \(a,z_1\),
where \(V=\{a=0\}\) and \(f=a^A z_1^B\epsilon\) with \(A>0\), \(B>0\)
and \(\epsilon\) a unit. Extend \(h\) from \(V\) to a regular local
function \(\widetilde h\) divisible by \(z_1\), which is possible
because \(h\) vanishes at the crossing. For the local one-form
\(\eta=\widetilde h\, da/a\), the logarithmic residue of \(d\eta\)
is \(-dh\) with the convention \(\operatorname{Res}_V(da/a\wedge\alpha)
=\alpha|_V\). Thus adding \(\nabla\eta\) to \(\omega\) removes the
ordinary simple pole along \(V\). The coefficient of the additional
\(\partial_b\)-term, up to sign, is
\[
 \widetilde h\,\frac{f_{z_1}}{a}
  =\widetilde h\,a^{A-1}z_1^{B-1}(B\epsilon+z_1\epsilon_{z_1}).
\]
It is a regular multiple of
\(f_a=a^{A-1}z_1^B(A\epsilon+a\epsilon_a)\).
At points of \(V\) away from the crossings take \(f=a^A\epsilon\);
the same divisibility follows directly without requiring a zero
of \(\widetilde h\). The graph relation
\[
 f_a\partial_b\delta_f=-\partial_a\delta_f
\]
puts the extra term in \(V^{>0}\), since \(\partial_a\) is a source
operator. At \(H\cap V\) take \(\eta=\widetilde h\,df/f\), extending
\(h\) as a meromorphic function with poles only along \(H\).
This gives the same residue, while \(df\wedge\eta=0\). The corrected
top coefficient then has only its allowed horizontal pole.

On open sets off \(V\) use \(\eta=0\). These local one-forms give a Čech
zero-cochain in de Rham degree one. With the total differential
\(d_{\mathrm{tot}}=d_{\check{C}}+(-1)^p\nabla\) on Čech degree \(p\),
\(\omega+d_{\mathrm{tot}}\eta\) represents the original graph class.
Its top components are \(\omega+\nabla\eta_i\) and belong to
\(V^{>0}\). The components in Čech degree one are
\(\eta_j-\eta_i\). These have zero residue on \(V\), hence have poles only
along the allowed horizontal divisors near the central fibre. They have
no \(\partial_b\) factor and lie in \(V^{>0}\) by the same local tests.
This cancellation uses the *same* primitive \(h\) on every chart.

Take a finite cover over the proper central fibre and a common positive
filtration index. The filtered relative de Rham pushforward theorem
puts this Čech class in \(V^{>0}\) after proper direct image
([Chen–Dirks–Mustaţă, Theorem 8.1](https://arxiv.org/pdf/2402.07816v1)).
The [proper-volume Fourier argument](../../../unpublished/primitive-balanced-image-spectrum/proper-volume.tex)
now excludes a nonzero polynomial exponential period. Čech correction
preserves the period because the original and corrected cocycles are
cohomologous in the localized graph complex, and the compact sphere
cycle avoids the boundary poles. Analytic specialization of the same
relative de Rham relations integrates these coboundaries to zero.
The constant term of the normalized sphere period is one, so its
positive coefficients cannot vanish eventually.

The retained degree-six off-center seed has
\(g(w)=e_4w^4+b_2w^2\) with \(b_2\ne0\) at this boundary. Its three
distinct finite roots place it outside the proposition. For several
roots, exactness of the residue by itself gives no common primitive
vanishing at every crossing, and the above filtered correction does
not resolve that case.
