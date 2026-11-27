# Three-variable Keller maps with a low-degree coordinate component

**Theorem ThreeVariableKellerBoundedCoordinateCompletion.** Let \(F\)
be a three-variable polynomial Keller map over a characteristic-zero
field \(k\), and let \(d\) be its ordinary total degree, the maximum
of its coordinate degrees. Suppose a nonconstant \(k\)-linear
combination \(R\) of its coordinates is a component of a polynomial
automorphism \(\Theta\) over \(k\). Write \(e=\deg(\Theta^{-1})\).
If \(de\le100\), then \(F\) is a polynomial automorphism.

After an invertible linear target change, write \(F=(P,Q,R)\).
Permute the coordinates of \(\Theta\) to put \(R\) last. Then
\[
 F\circ\Theta^{-1}=(\widetilde P,\widetilde Q,t)
 \qquad\text{and}\qquad
 \max(\deg\widetilde P,\deg\widetilde Q)\le de.
\]
This map still has nonzero constant Jacobian determinant, by the chain
rule and invertibility of \(\Theta\). Its restriction to every plane
\(t=t_0\) has the same nonzero constant two-variable Jacobian
determinant. Specialization does not increase the degree.

First work over \(\mathbb C\). On each such plane the
[retained Moh degree-100 theorem](../../BREAKTHROUGH_TARGET_ASSESSMENT.md#four-variable-hessian-conjecture-and-planar-degree-bounds)
gives a polynomial automorphism. The third coordinate \(t\) separates
the planes, so \(F\circ\Theta^{-1}\) is injective on complex points.
[Ax–Grothendieck](https://msp.org/pjm/1969/31-1/pjm-v31-n1-p01-p.pdf),
Corollary 1, gives surjectivity. The Keller determinant makes the map
étale. Its diagonal in the fibre product is both an open and a closed
immersion, and the complement has no complex point by injectivity.
A nonempty finite-type complex scheme has a closed complex point, so
the diagonal is an isomorphism. The map is therefore universally
injective. A surjective étale universally injective map is an
isomorphism by the
[open immersion criterion](https://stacks.math.columbia.edu/tag/025G).
This proves invertibility of \(F\).

For general \(k\), take a finitely generated subfield \(k_0\subset k\)
containing the coefficients of \(F,\Theta,\Theta^{-1}\) and of the
linear target change. Embed \(k_0\) in \(\mathbb C\). After subtracting
the value of \(F\) at the origin, its unique formal inverse at the
origin has coefficients in \(k_0\), by recursive inversion of its
invertible linear term. Over \(\mathbb C\) it is the polynomial inverse
just proved, so the formal series over \(k_0\) terminates. The polynomial
inverse extends to \(k\) and to the original target translation.

**Corollary DegreeSixKellerQuadraticComponent.** A three-variable
polynomial Keller map over a characteristic-zero field of ordinary
total degree at most six is invertible if a nonconstant linear
combination \(R\) of its coordinates has degree at most two. The same
conclusion holds if such a combination has an explicitly supplied
polynomial coordinate completion with inverse degree at most five.

For the quadratic assertion, write \(R=q+\ell+c\), where \(q\) is
homogeneous quadratic and \(\ell\) is linear. The differential of
\(R\) can never vanish, since it is a nonzero constant combination
of rows of an invertible Jacobian matrix. The linear map
\(v\mapsto \mathrm dq(v)\) has image the annihilator of the radical
of the symmetric form associated to \(q\). Consequently \(\ell\)
does not vanish on that radical, or \(\mathrm dq(v)+\ell=0\)
would have a solution. Choose a radical direction on which \(\ell\)
is one and call its coordinate \(z\). In a complementary linear
frame, \(R=z+s(x,y)\) with \(\deg s\le2\). The triangular completion
\((x,y,R)\) has inverse degree at most two. The theorem applies with
\(de\le12\). For the supplied inverse-degree-five completion it applies
with \(de\le30\).

This is the bounded-degree planar reduction used by Wang and Yan in
[arXiv:2610.36495v1](https://arxiv.org/pdf/2610.36495v1), Lemma 3.12.
The explicit completion hypothesis gives the conclusion irrespective
of the degree of \(R\) within the total-degree bound. The cubic
classification in their Proposition 3.16 has separate quartic
hypotheses.
