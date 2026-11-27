# The three-stage embedding at conductor fourteen

**Theorem ConductorFourteenThreeStageExclusion.** Let \((P,Q)\)
be a complex polynomial plane Keller map in a source frame with
\(P\) monic of positive longitudinal degree. Its generic coordinate
order cannot have reduced characteristic sequence \((12,8,6,3)\).
The degree-six and degree-three approximate roots together force
the coefficient below the leading term of the degree-eight generator
to vanish in the frame centred at the degree-three element. This
removes the constant term of the parameter Jacobian.

## The two approximate-root identities

Use the [degree-fourteen classification](CONDUCTOR_FOURTEEN_CHARACTERISTIC_SEQUENCES.md)
and its [approximate-root input](sources/assi-garcia-sanchez-1407.0490v1.pdf),
namely Assi–García-Sánchez Proposition 2 and the setup for its approximate
roots. Over \(K=\mathbb C(y)\) choose monic reduced generators \(V,W\)
of degrees \(12,8\) in the longitudinal variable. Let \(S,T\) be the
successive approximate roots of longitudinal degrees six and three.
As polynomials in \(W\) over \(K[V]\), they are monic of degrees
three and six, respectively.

The degree-six identity has \(S=W^3+\sum_{i<3} a_i(V)W^i\).
The degrees of the latter terms are distinct modulo \(12\) and
only \(W^3\) and a multiple of \(V^2\) can have degree \(24\).
It follows that \(a_0\) is quadratic, \(a_1\) affine and \(a_2\)
constant. Complete the square in \(V\) by a linear polynomial
in \(W\), then translate \(W\) to depress its cube. With monic
generators still denoted \(V,W\), there is an identity
\[
                        V^2=W^3+kW+R(T),
 \qquad k\in K,\quad R(T)\in K[T],\quad \deg_T R=2 .
\tag{1}
\]
Here \(S\) is quadratic in \(T\). To see this, the characteristic
sequence generates the semigroup \(\langle3,8\rangle\), and its
only degrees below eight are \(0,3,6\). For any element of degree
six, subtract multiples of \(T^2,T,1\) in descending order of
degree. This expresses \(S\) as a quadratic with nonzero
leading coefficient.

For the degree-three approximate root, divide the monic degree-six
polynomial in \(W\) by the monic degree-three \(S\) over \(K[V]\).
The result has the form \(S^2+\sum_{i<3,j<2} b_{ij}(V)W^iS^j\).
The longitudinal degrees \(8i+6j\) have distinct residues modulo
\(12\) for these six choices of \((i,j)\). Thus no summand of
degree greater than \(12\) can cancel. The only terms of degree at
most \(12\) are scalar multiples of \(1,V,W,S\), and \(V\) has
nonzero coefficient because it must cancel the highest term of
\(S^2\). This gives the additional relation
\[
                        V=A(T)+B W ,
 \qquad A\in K[T],\quad \deg A=4,\quad B\in K .
\tag{2}
\]
The square completion in \(V\) and translation in \(W\) preserve
this form. We can center the longitudinal variable \(u\) and
affinely normalize \(T\) so \(T=u^3+a u\). Equations (1) and (2)
are preserved under this change.

## The coefficient below the leading term of \(W\)

The field \(K(T,W)\) equals \(K(u)\) by the coprime polynomial
degrees \(3,8\). Over \(K(T)\) the element \(W\) therefore has
degree three, and (1) and (2) give its monic relation
\[
             W^3-B^2W^2+(k-2AB)W+R-A^2=0.
\]
In particular \(\operatorname{Tr}(W)=B^2\) and
\(\operatorname{Tr}(W^2)=B^4-2(k-2AB)\) has \(T\)-degree at
most four.

Use the polynomial basis \(1,u,u^2+2a/3\) over \(K[T]\) and write
\[
      W=w_0(T)+p(T)u+q(T)(u^2+2a/3).
\]
The three longitudinal degree residues modulo three are distinct,
so \(\deg p\le2\) and \(q\) is monic of degree two. The two
nonconstant basis elements have trace zero, since the roots of
\(u^3+a u-T\) have power traces
\(\operatorname{Tr}(u)=0\), \(\operatorname{Tr}(u^2)=-2a\),
\(\operatorname{Tr}(u^3)=3T\) and
\(\operatorname{Tr}(u^4)=2a^2\). Hence \(w_0=B^2/3\).
The mixed term of \(\operatorname{Tr}(W^2)\) is \(6T p(T)q(T)\),
whereas the trace of each separate square has degree at most four.
The coefficient of \(T^5\) forces the quadratic coefficient of
\(p(T)\) to vanish. Since \(T=u^3+a u\), the only possible
contribution to the coefficient of \(u^7\) in \(W\) was that
quadratic coefficient. Thus \(W=u^8+O(u^6)\) at infinity.

## The parameter Jacobian

Apply the
[determinant-one normalization](CONDUCTOR_TWELVE_CABLE_EXCLUSIONS.md#1-common-source-and-target-identities).
If the determinant of the target change is \(\lambda\in K^\times\),
take \(u=\rho z\) with \(\rho^{20}=\lambda\) and scale \(W,V\)
by \(\rho^{-8},\rho^{-12}\), respectively. Normalize \(T\) by
\(\rho^{-3}\). These operations preserve (1), (2), and both
coefficient vanishings. We retain the names \(W,V,T,R\) for the
normalized expressions. In particular \(W=z^8+O(z^6)\) and
\(R=r_2z^6+O(z^4)\), since \(R\) is quadratic in the depressed
\(T=z^3+a z\). In the algebraic parameter field their parameter
Jacobian satisfies
\[
 W_zV_y-W_yV_z
      = J\rho+\partial_z H(W,V),
 \qquad J\in\mathbb C^\times ,
\tag{3}
\]
where \(H\) is polynomial in the two target variables. The
parameter derivative of the determinant-one target change is
a divergence-free polynomial field. Its Hamiltonian gives the
last term of (3).

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

The semigroup \(\langle3,8\rangle\) implies that the evaluated
\(H(W,V)\) in (3), of degree at most four, belongs to the span
of \(1,T\). Comparing the degree-three coefficients gives \(k'=0\).
Comparing degree two then makes the coefficient of \(T\) in \(H\)
zero. The constant term of the parameter Jacobian is now zero,
in contradiction to \(J\rho\ne0\) in (3).
