# Three-stage embeddings with first degrees twelve, eight and six

**Theorem KellerTwelveEightSixThreeStageExclusion.** Let
\((P,Q)\) be a complex polynomial plane Keller map in a source
frame with \(P\) monic of positive longitudinal degree. Its
generic coordinate order cannot have reduced characteristic
sequence \((12,8,6,r)\) for \(r=3,7,9\). The successive
approximate roots force two generator centroids to coincide.
Centering the longitudinal coordinate at that point removes
the constant term of the parameter Jacobian.

## The two approximate-root identities

Use the [degree-fourteen classification](CONDUCTOR_FOURTEEN_CHARACTERISTIC_SEQUENCES.md)
and its [approximate-root input](sources/assi-garcia-sanchez-1407.0490v1.pdf),
namely Assi–García-Sánchez Proposition 2 and the setup for its approximate
roots. Over \(K=\mathbb C(y)\) choose monic reduced generators \(V,W\)
of degrees \(12,8\) in the longitudinal variable. Let \(S,T\) be the
successive approximate roots of longitudinal degrees six and \(r\).
As polynomials in \(W\) over \(K[V]\), they are monic of degrees
three and six, respectively.

The degree-six identity has \(S=W^3+\sum_{i<3} a_i(V)W^i\).
The degrees of the latter terms are distinct modulo \(12\) and
only \(W^3\) and a multiple of \(V^2\) can have degree \(24\).
It follows that \(a_0\) is quadratic, \(a_1\) affine and \(a_2\)
constant. Complete the square in \(V\) by a linear polynomial
in \(W\), then translate \(W\) to depress its cube. Rescale
\(S\) to a monic polynomial \(S_0\), and keep the names \(V,W\)
for the monic target generators of degrees \(12,8\). There is
an identity
\[
                       V^2=W^3+kW+\eta S_0+C .
\tag{1}
\]

For the last approximate root, divide the monic degree-six
polynomial \(T\) in \(W\) by the monic degree-three \(S\) over
\(K[V]\) in the original target coordinates. The result is
\(S^2+\sum_{i<3,j<2} b_{ij}(V)W^i S^j\). The longitudinal
degrees \(8i+6j\) have distinct residues modulo \(12\) for these
six choices of \((i,j)\). No summand of degree greater than
\(12\) can cancel. The only terms of degree at most \(12\)
are scalar multiples of \(1,V,W,S\), and \(V\) has nonzero
coefficient because it must cancel the highest term of \(S^2\).
In the new coordinates of (1) the same expression for \(T\)
therefore lies in the span of \(S_0^2,V,W,S_0,1\). Since
\(\deg T\le9\), it follows by monicity that
\[
                       V=S_0^2+O(u^9)
\tag{2}
\]
at longitudinal infinity.

## The common centroid and the parameter Jacobian

Write \(w_7,v_{11},s_5\) for the coefficients below the leading
terms of \(W,V,S_0\), respectively. From (1) at degree \(23\)
one gets \(2v_{11}=3w_7\), and from (2) one gets
\(v_{11}=2s_5\). Thus \(w_7/8=s_5/6\).
Center the source by translating \(x\) by \(h(y)=-w_7/8\).
This removes both the degree-seven term of \(W\) and the
degree-five term of \(S_0\).

Now apply the
[determinant-one normalization](CONDUCTOR_TWELVE_CABLE_EXCLUSIONS.md#1-common-source-and-target-identities).
If the determinant of the target change to \((W,V)\) is
\(\lambda\in K^\times\), take \(x=h(y)+\rho z\) with
\(\rho^{20}=\lambda\) and scale \(W,V,S_0\) by
\(\rho^{-8},\rho^{-12},\rho^{-6}\), respectively. Equation (1)
keeps its form with rescaled coefficients and both coefficient
vanishings persist. Put \(R=\eta S_0+C\) in this frame, so
\(W=z^8+O(z^6)\) and \(R=r_2 z^6+O(z^4)\). The parameter
Jacobian satisfies
\[
 W_zV_y-W_yV_z
      =J\rho+\partial_z H(W,V),
 \qquad J\in\mathbb C^\times .
\tag{3}
\]
Here \(H\) is polynomial in the target variables. The parameter
derivative of the determinant-one target change is a
divergence-free polynomial field, and its Hamiltonian gives
the last term of (3).

Expand (1) in the monic Laurent branch at \(z=\infty\).
\[
 V=W^{3/2}+\tfrac12 k W^{-1/2}
               +\tfrac12 R W^{-3/2}+O(z^{-20}).
\]
The first term contributes zero to the left side of (3). The
second contributes \(k'(\sqrt W)_z\), whose polynomial part has
degree at most three and top coefficient \(4k'\). The third
contributes \(\tfrac12 W^{-3/2}(W_zR_y-W_yR_z)\).
Here \(W_zR_y-W_yR_z=8r_2' z^{13}+O(z^{11})\) and
\(W^{-3/2}=z^{-12}(1+O(z^{-2}))\). Its polynomial part is
\(4r_2'z\), with zero constant. The last Laurent tail only
contributes negative powers. All comparisons hold with the
extended derivation on the finite coefficient field.

The pole semigroup is \(\langle3,8\rangle\) for \(r=3\)
and \(\langle6,8,r\rangle\) for \(r=7,9\). Thus the evaluated
\(H(W,V)\) in (3), of degree at most four, belongs to the span
of \(1,T\) for \(r=3\) and is constant for \(r=7,9\).
Comparing degree three gives \(k'=0\). Comparing degree two
eliminates the possible \(T\) coefficient. The constant term
in (3) is now zero on the left and \(J\rho\ne0\) on the right.
