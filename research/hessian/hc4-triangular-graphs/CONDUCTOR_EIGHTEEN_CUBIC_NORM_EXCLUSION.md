# The cubic norm excludes the conductor-eighteen sequence (10,6,3)

**Theorem ConductorEighteenCubicOrderTenSixThreeExclusion.** Over a
characteristic-zero field \(K\), a polynomial order \(K[W,V]\subset K[z]\)
with fraction field \(K(z)\) and normalization \(K[z]\) cannot have
reduced one-place characteristic sequence \((10,6,3)\).
In particular this sequence cannot occur for the exact directional
order of a polynomial plane Keller map.

The contradiction uses only the norm over the degree-three element
of the order. It applies before using a parameter equation.

## The trace-zero equation

Pass to an algebraic closure of \(K\) if necessary. Degree-filtered
dimensions and the fraction-field normalization persist under this
extension. Let the reduced generators \(W,V\) be monic of degrees six
and ten. Their pole semigroup is \(\langle3,10\rangle\), and the
order has a degree-three element. Normalize one such element by
constant scaling and translation of the source to a monic depressed
cubic. The filtered order through degree six has basis consisting of
one, this cubic, and its square. Expressing \(W\) in that basis and
completing the square, with a constant target translation, gives
\[
                       W=U^2,\qquad U=z^3+az+h .
\]
The degree-three approximate root remains affine in \(U\), with
nonzero slope. The approximate-root input is the one-place theorem
of Assi–García-Sánchez, Propositions 2 and 4, used with the coefficient
degree bounds and normalization from the
[previous cubic-order proof](CONDUCTOR_FOURTEEN_CUBIC_ORDER_ALTERNATE_EMBEDDINGS.md).
[Retained primary source](sources/assi-garcia-sanchez-1407.0490v1.pdf)

In this case the approximate root is monic of degree five in \(W\)
and has source degree three. For weights \((6,10)\), the exponents
\(0,1,2,3,4\) of \(W\) give distinct weights modulo ten. Thus terms
above weight thirty cannot cancel, and the top terms are \(W^5-V^3\).
The coefficient of \(V^2\) is affine in \(W\). A triangular target
change adding an affine function of \(W\) to \(V\) depresses the cubic
and preserves its monic leading term. The equation of the actual
normalized pair takes the form
\[
 V^3+A(U)V+B(U)=0,\qquad A(-U)=A(U),\qquad
             B(U)-B(-U)=2\gamma U,\quad \gamma\ne0 .        \tag{1}
\]
Indeed the target coefficients are polynomials in \(W=U^2\),
and the approximate root itself is affine in \(U\).

The extension \(K(z)/K(U)\) has degree three. The element \(V\)
of source degree ten does not belong to \(K(U)\), since the pole
order of a rational function of \(U\) at infinity is divisible
by three. Hence (1) is its minimal polynomial, and \(V\) has
trace zero. In the \(K[U]\)-basis \(1,z,z^2\),
\(\operatorname{Tr}(z)=0\) and \(\operatorname{Tr}(z^2)=-2a\).
It follows that
\[
 \begin{aligned}
 V&=q(U)z+p(U)(z^2+2a/3),\\
 q(U)&=U^3+q_2U^2+q_1U+q_0,\qquad
 p(U)=p_2U^2+p_1U+p_0 .
 \end{aligned}
\]
These degree bounds and the leading coefficient of \(q\) follow
from the three distinct degree residues modulo three in the basis.

## Coefficients of the norm

The norm formulas in the
[conductor-fourteen trace-zero calculation](CONDUCTOR_FOURTEEN_CUBIC_ORDER_ALTERNATE_EMBEDDINGS.md#a-sextic-lower-generator-and-an-octic-upper-generator)
apply to these \(p,q\):
\[
\begin{aligned}
 A &=a q^2-3(U-h)pq-a^2p^2/3,\\
 B &=-(U-h)^2p^3+(U-h)ap^2q-(U-h)q^3
          -2a^3p^3/27-2a^2pq^2/3.
\end{aligned}                                             \tag{2}
\]
Write \(A_i=[U^i]A\), \(B_i=[U^i]B\). The successive equations
\(B_9=A_5=B_7=0\) give
\[
 q_2=h/3,\qquad
 p_1=2ah/9+2hp_2/3,\qquad
 q_0=h(-4a^3-12a^2p_2+8h^2+27q_1)/81.
\]
If \(h=0\), then \(q\) is odd and \(p\) even, so \(B\) is even
by (2). This contradicts \(\gamma\ne0\). Thus \(h\ne0\).
Now \(A_3=0\) gives
\[
 p_0=(4a^4+12a^3p_2-17ah^2-27a q_1-15h^2p_2)/81 .
\]
Set \(z=\eta w\) with \(\eta^3=h\), and scale the source polynomials
\(U,W,V\) by \(\eta^{-3},\eta^{-6},\eta^{-10}\), respectively.
This takes \(h\) to one and preserves monicity and (1), with
rescaled coefficients. Use the original symbols for these
new coefficients.

After the displayed comparisons, \(B_5=0\) gives
\(q_1=13a^3/27+a^2p_2-5/9\). Substitution gives
\[
 \begin{split}
 A_1&=32(5a+6p_2)/2187,\\
 B_3&=32(-15a^3-28a^2p_2+3ap_2^2+18p_2^3+4)/6561 .
 \end{split}
\]
Thus \(A_1=0\) forces \(p_2=-5a/6\), where
\(B_3=128/6561\ne0\). The contradiction with (1) proves the theorem.

The [exact checker](verify_conductor_eighteen_cubic_norm_exclusion.py)
replays the norm and all these coefficient eliminations.
